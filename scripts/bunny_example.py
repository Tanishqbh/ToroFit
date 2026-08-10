"""Run the merged 05_bunny notebook workflow and export SVG figures."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from dreimac import ToroidalCoords
from matplotlib.text import Text
from scipy.special import softmax

import torofit as TF
from torofit.partition_utils import total_variation_scores

try:
    import cmasher as cmr  # noqa: F401  Registers cmr.* colormaps with Matplotlib.
except ImportError:
    cmr = None


COLOR_MAP_1 = "hsv"
COLOR_MAP_2 = "cmr.seasons_s"
FALLBACK_COLOR_MAP_2 = "twilight_shifted"
REPO_ROOT = Path(__file__).resolve().parents[1]


def get_color_maps() -> list:
    color_map_2 = COLOR_MAP_2 if COLOR_MAP_2 in plt.colormaps() else FALLBACK_COLOR_MAP_2
    return [plt.colormaps[COLOR_MAP_1], plt.colormaps[color_map_2]]


def generate_s1_s2_s1_point_cloud(
    n_points: int = 1000,
    n_circle: int = 500,
    sphere_radius: float = 8.0,
    circle1_radius: float = 5.0,
    circle2_radius: float = 5.1,
    wedge_point: np.ndarray | None = None,
    circle1_normal: np.ndarray | None = None,
    circle2_normal: np.ndarray | None = None,
    noise: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate the notebook's S^1 wedge S^2 wedge S^1 bunny point cloud."""
    phi = np.random.uniform(0, 2 * np.pi, n_points)
    z = np.random.uniform(-1, 1, n_points)
    r = sphere_radius * np.sqrt(1 - z**2)
    sphere_points = np.column_stack(
        (r * np.cos(phi), r * np.sin(phi), sphere_radius * z)
    )

    if wedge_point is None:
        wedge_point = np.array([0, 0, sphere_radius])
    if circle1_normal is None:
        circle1_normal = np.array([0, 1, 0])
    if circle2_normal is None:
        circle2_normal = np.array([np.sin(np.pi / 4), np.cos(np.pi / 4), 0])

    angles = np.linspace(0, 2 * np.pi, n_circle, endpoint=False)

    def make_circle(p0, normal, radius, angles):
        p0 = np.array(p0, dtype=float)
        n = np.array(normal, dtype=float)
        n = n / np.linalg.norm(n)
        aux = np.array([1, 0, 0]) if abs(np.dot(n, [1, 0, 0])) < 0.9 else np.array([0, 1, 0])
        e1 = np.cross(n, aux)
        e1 = e1 / np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        e2 = e2 / np.linalg.norm(e2)
        center = p0 + radius * e1
        return np.array(
            [center + radius * (np.cos(t) * e1 + np.sin(t) * e2) for t in angles]
        )

    circle1 = make_circle(wedge_point, circle1_normal, circle1_radius, angles)
    circle2 = make_circle(wedge_point, circle2_normal, circle2_radius, angles)
    point_cloud = np.vstack([sphere_points, circle1, circle2])
    point_cloud += np.random.normal(0, noise, point_cloud.shape)
    return point_cloud, sphere_points, circle1, circle2


def partition_points(
    tc: ToroidalCoords,
    theta: np.ndarray,
    *,
    knn: int = 20,
    thresholds: tuple[float, float] = (0.505, 0.505),
    low_score_cutoff: float = 0.02,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Partition data points and compute notebook-style softmax weights."""
    part1 = total_variation_scores(0, tc, theta, K=knn)
    part2 = total_variation_scores(1, tc, theta, K=knn)

    raw_scores = np.column_stack([part1, part2])
    softmax_scores = softmax(raw_scores, axis=1)
    filter_for_low_raw_scores = np.linalg.norm(raw_scores, axis=1) > low_score_cutoff
    assigned_part = -np.ones(raw_scores.shape[0], dtype=int)

    for i, threshold in enumerate(thresholds):
        confident_idx = np.where(
            (softmax_scores[:, i] >= threshold) & filter_for_low_raw_scores
        )[0]
        assigned_part[confident_idx] = i

    parts = []
    weights = []
    landmarks = tc._idx_land
    for i in range(raw_scores.shape[1]):
        good_idx = np.where(assigned_part == i)[0]
        print(f"Partition {i + 1}, threshold: {thresholds[i]}")
        print(f"  confident landmarks in partition: {len(good_idx)}")

        if len(good_idx) == 0:
            parts.append(np.array([], dtype=int))
            weights.append(np.array([], dtype=float))
            continue

        _ = landmarks[good_idx]
        mems_for_good = tc._dist_land_data[good_idx, :]
        in_cover = mems_for_good < tc._r_cover
        part_pts = np.where(np.any(in_cover, axis=0))[0]
        parts.append(part_pts)

        scores_i = softmax_scores[good_idx, i]
        masked = np.where(in_cover[:, part_pts], scores_i[:, None], 0.0)
        weights.append(masked.max(axis=0))

    print(f"Ignored low-score landmarks: {np.sum(assigned_part == -1)}")
    print(f"Partition 1 data points: {len(parts[0])}")
    print(f"Partition 2 data points: {len(parts[1])}")
    return parts, weights


def assembly(
    fit,
    theta,
    K,
    wedge_point_coords,
    samples=100,
    wedge_point=np.array([0, 0, 8]),
):
    if theta is None:
        theta = np.linspace(0, 2 * np.pi, num=samples).reshape(-1, 1)
    image = TF.Toroidal_Evaluate(fit, theta + wedge_point_coords, K).real
    offset = wedge_point - TF.Toroidal_Evaluate(
        fit, wedge_point_coords.reshape(1, -1), K=K
    ).real
    image += offset
    return image


class ColorMapHandler:
    """Legend handler that draws a horizontal colormap bar."""

    def __init__(self, cmap_name, steps=100, bar_height_frac=0.7, y_shift_frac=0.12):
        if cmap_name not in plt.colormaps():
            cmap_name = FALLBACK_COLOR_MAP_2
        self.cmap = plt.colormaps[cmap_name]
        self.steps = steps
        self.bar_height_frac = bar_height_frac
        self.y_shift_frac = y_shift_frac

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        bar_height = self.bar_height_frac * height
        y_bar = y0 + 0.5 * (height - bar_height) - self.y_shift_frac * height
        dx = width / self.steps
        first_rect = None
        for i in range(self.steps):
            rect = patches.Rectangle(
                (x0 + i * dx, y_bar),
                dx,
                bar_height,
                facecolor=self.cmap(i / self.steps),
                edgecolor="none",
            )
            if first_rect is None:
                first_rect = rect
            handlebox.add_artist(rect)
        return first_rect


class TextOnlyHandler:
    """Legend handler that centers custom text."""

    def __init__(self, label, color="black", y_shift_frac=0.0):
        self.label = label
        self.color = color
        self.y_shift_frac = y_shift_frac

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        text = Text(
            x=x0 + width / 2,
            y=y0 + height / 2 - self.y_shift_frac * height,
            text=self.label,
            ha="center",
            va="center",
            color=self.color,
            fontsize=fontsize,
        )
        handlebox.add_artist(text)
        return text


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    max_range = 0.5 * max(
        abs(x_limits[1] - x_limits[0]),
        abs(y_limits[1] - y_limits[0]),
        abs(z_limits[1] - z_limits[0]),
    )
    x_mid = np.mean(x_limits)
    y_mid = np.mean(y_limits)
    z_mid = np.mean(z_limits)
    ax.set_xlim3d([x_mid - max_range, x_mid + max_range])
    ax.set_ylim3d([y_mid - max_range, y_mid + max_range])
    ax.set_zlim3d([z_mid - max_range, z_mid + max_range])


def clean_3d_axes(ax):
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("#dddddd")
    ax.yaxis.pane.set_edgecolor("#dddddd")
    ax.zaxis.pane.set_edgecolor("#dddddd")


def plot_theta_colored_line(ax, curve, colors, linewidth=2.2, alpha=0.95):
    for j in range(curve.shape[0] - 1):
        ax.plot(
            curve[j : j + 2, 0],
            curve[j : j + 2, 1],
            curve[j : j + 2, 2],
            color=colors[j],
            linewidth=linewidth,
            alpha=alpha,
        )


def configure_matplotlib():
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 500,
            "savefig.facecolor": "white",
            "font.family": "DejaVu Serif",
            "font.size": 10.5,
            "axes.labelsize": 11.5,
            "axes.titlesize": 12.5,
            "legend.fontsize": 10,
        }
    )


def save_fit_overlay(output_dir, point_cloud, curve1, curve2, theta_plot, color_maps):
    text_proxy_1 = patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none")
    text_proxy_2 = patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none")
    proxy_handle1 = patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none")
    proxy_handle2 = patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none")

    fig = plt.figure(figsize=(12.2, 5.6), facecolor="white")
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    theta_colors = [
        color_maps[0](theta_plot[:, 0] / (2 * np.pi)),
        color_maps[1](theta_plot[:, 0] / (2 * np.pi)),
    ]

    for ax in [ax1, ax2]:
        ax.set_facecolor("white")
        ax.scatter(
            point_cloud[:, 0],
            point_cloud[:, 1],
            point_cloud[:, 2],
            s=5,
            c="gray",
            alpha=0.2,
            linewidth=0,
            rasterized=True,
        )
        plot_theta_colored_line(ax, curve2, theta_colors[1], linewidth=3, alpha=0.95)
        plot_theta_colored_line(ax, curve1, theta_colors[0], linewidth=3, alpha=0.95)
        ax.set_xlabel("$x$", labelpad=2)
        ax.set_ylabel("$y$", labelpad=2)
        ax.set_zlabel("$z$", labelpad=2)
        set_axes_equal(ax)
        clean_3d_axes(ax)

    ax1.view_init(elev=16, azim=-68)
    ax2.view_init(elev=29, azim=22)
    ax1.grid(True, color="#dddddd", linewidth=0.5, alpha=0.8)
    ax2.grid(True, color="#dddddd", linewidth=0.5, alpha=0.8)
    fig.legend(
        [text_proxy_1, text_proxy_2, proxy_handle1, proxy_handle2],
        ["", "", "ToroFit 1", "ToroFit 2"],
        handler_map={
            text_proxy_1: TextOnlyHandler("Torofit output", y_shift_frac=0.55),
            text_proxy_2: TextOnlyHandler(r"$F: S^1 \vee S^1 \to X^{(\alpha)}$", y_shift_frac=0.55),
            proxy_handle1: ColorMapHandler(COLOR_MAP_1, bar_height_frac=0.6, y_shift_frac=0.55),
            proxy_handle2: ColorMapHandler(COLOR_MAP_2, bar_height_frac=0.6, y_shift_frac=0.55),
        },
        loc="lower center",
        bbox_to_anchor=(0.5, 0.03),
        ncol=2,
        columnspacing=2.2,
        handlelength=6,
        handleheight=2,
        markerfirst=False,
        frameon=True,
        alignment="center",
    )
    fig.subplots_adjust(top=0.88, bottom=0.12, wspace=0.04)
    fig.savefig(output_dir / "05_bunny_fits.svg", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def save_point_cloud_with_persistence(output_dir, point_cloud, dgms):
    from persim import plot_diagrams

    fig = plt.figure(figsize=(10.6, 5.6), facecolor="white")
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1], wspace=0.05, hspace=0.05)
    ax = fig.add_subplot(gs[0], projection="3d")
    diag_ax = fig.add_subplot(gs[1])
    ax.scatter(
        point_cloud[:, 0],
        point_cloud[:, 1],
        point_cloud[:, 2],
        s=3.5,
        c="gray",
        alpha=0.65,
        linewidth=0,
        rasterized=True,
    )
    ax.set_xlabel("$x$", labelpad=2)
    ax.set_ylabel("$y$", labelpad=2)
    ax.set_zlabel("$z$", labelpad=2)
    set_axes_equal(ax)
    clean_3d_axes(ax)
    ax.view_init(elev=30, azim=-10)
    ax.grid(True, color="#dddddd", linewidth=0.5, alpha=0.8)

    plot_diagrams(dgms, ax=diag_ax, size=30)
    diag_ax.legend(loc="lower right", frameon=True, fontsize=9)
    if len(dgms) > 1 and len(dgms[1]) >= 2:
        perm = np.argsort(dgms[1][:, 1] - dgms[1][:, 0])
        target_class = dgms[1][perm[-2:]]
        axins = diag_ax.inset_axes([0.7, 0.5, 0.25, 0.25])
        plot_diagrams(dgms, ax=axins, size=30, legend=False)
        axins.set_xticks([])
        axins.set_yticks([])
        x_target, y_target = target_class[0, 0], target_class[0, 1]
        margin = 0.5
        axins.set_xlim(x_target - margin, x_target + margin)
        axins.set_ylim(y_target - margin, y_target + margin)
        diag_ax.indicate_inset_zoom(axins, edgecolor="black")

    fig.tight_layout()
    fig.savefig(output_dir / "05_bunny_point_cloud.svg", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def save_toroidal_coordinates(output_dir, point_cloud, theta, color_maps):
    proxy_handle1 = patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none")
    proxy_handle2 = patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none")
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8), facecolor="white")
    for i, ax in enumerate(axes):
        ax.scatter(
            point_cloud[:, 0],
            point_cloud[:, 1],
            c=color_maps[i](theta[i] / (2 * np.pi)),
            s=10,
            alpha=0.75,
            linewidth=0,
            rasterized=True,
        )
        ax.set_aspect("equal")
        ax.axis("off")

    fig.tight_layout()
    fig.legend(
        [proxy_handle1, proxy_handle2],
        ["Coord 1 (HSV Colormap)", "Coord 2 (Seasons Colormap)"],
        handler_map={
            proxy_handle1: ColorMapHandler(COLOR_MAP_1),
            proxy_handle2: ColorMapHandler(COLOR_MAP_2),
        },
        loc="lower center",
        bbox_to_anchor=(0.5, 0.03),
        ncol=2,
        handlelength=4,
        handleheight=1.2,
    )
    fig.savefig(output_dir / "05_bunny_toroidal_coords.svg", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def save_cup_product_table(output_dir, point_cloud, landmarks):
    from torofit.plotting_utils import plot_cup_product_table
    from torofit.rips_complex_utils import get_nonzero_cup_products_dim_1

    tc = ToroidalCoords(point_cloud, n_landmarks=landmarks)
    _ = tc.get_coordinates(cocycle_idxs=[0, 1], perc=0.1)
    nonzero_cup_products = get_nonzero_cup_products_dim_1(
        tc, harmonic_class_idx=[0, 1], verbose=True
    )
    fig, ax_table = plot_cup_product_table(nonzero_cup_products, harmonic_class_idx=[0, 1])
    ax_table.set_title("")
    fig.savefig(output_dir / "05_bunny_cup_product_table.svg", bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)


def run(args):
    configure_matplotlib()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = REPO_ROOT / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    np.random.seed(args.seed)

    wedge_point = np.array([0, 0, 8])
    point_cloud, _, _, _ = generate_s1_s2_s1_point_cloud(
        n_points=args.sphere_points,
        n_circle=args.circle_points,
        wedge_point=wedge_point,
        circle1_normal=np.array([0, -1, 1]),
        circle2_normal=np.array([0, -1, -1]),
        noise=args.noise,
    )
    color_maps = get_color_maps()
    landmarks = len(point_cloud) if args.landmarks is None else args.landmarks

    print("Computing toroidal coordinates...")
    tc_bunny = ToroidalCoords(point_cloud, n_landmarks=landmarks)
    tc_total_theta = tc_bunny.get_coordinates(cocycle_idxs=[0, 1])

    parts, weights = partition_points(
        tc_bunny,
        tc_total_theta,
        knn=args.partition_knn,
        thresholds=(args.threshold_1, args.threshold_2),
        low_score_cutoff=args.low_score_cutoff,
    )

    print("Fitting weighted Fourier loops...")
    circle_1_fit = TF.fit(
        point_cloud[parts[0], :],
        tc_total_theta[0][parts[0]],
        K=args.bandwidth,
        lam=args.regularization,
        weights=weights[0],
    )
    circle_2_fit = TF.fit(
        point_cloud[parts[1], :],
        tc_total_theta[1][parts[1]],
        K=args.bandwidth,
        lam=args.regularization,
        weights=weights[1],
    )

    wedge_idx = np.argmin(np.linalg.norm(point_cloud - wedge_point, axis=1))
    wedge_point_coords = tc_total_theta[:, wedge_idx]
    theta_plot = np.linspace(0, 2 * np.pi, num=args.curve_samples).reshape(-1, 1)
    curve1 = assembly(circle_1_fit, theta_plot, K=args.bandwidth, wedge_point_coords=wedge_point_coords[0])
    curve2 = assembly(circle_2_fit, theta_plot, K=args.bandwidth, wedge_point_coords=wedge_point_coords[1])
    curve2[:, 0] = curve2[:, 0] - 0.05

    save_fit_overlay(output_dir, point_cloud, curve1, curve2, theta_plot, color_maps)
    save_toroidal_coordinates(output_dir, point_cloud, tc_total_theta, color_maps)

    if not args.skip_persistence:
        from ripser import ripser

        print("Computing persistence diagram...")
        dgms = ripser(point_cloud, maxdim=2, n_perm=args.persistence_landmarks)["dgms"]
        save_point_cloud_with_persistence(output_dir, point_cloud, dgms)

    if not args.skip_cup_product:
        print("Computing cup-product table...")
        save_cup_product_table(output_dir, point_cloud, max(1, len(point_cloud) // 3))

    print(f"SVG exports written to {output_dir}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="presentation_plots", help="Directory for SVG exports.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed matching the notebook.")
    parser.add_argument("--sphere-points", type=int, default=1000, help="Number of sphere samples.")
    parser.add_argument("--circle-points", type=int, default=500, help="Number of samples per circle.")
    parser.add_argument("--landmarks", type=int, default=None, help="DREiMac landmarks. Defaults to len(point_cloud), matching the notebook.")
    parser.add_argument("--partition-knn", type=int, default=20, help="KNN value for total-variation partition scores.")
    parser.add_argument("--threshold-1", type=float, default=0.505, help="Softmax threshold for partition 1.")
    parser.add_argument("--threshold-2", type=float, default=0.505, help="Softmax threshold for partition 2.")
    parser.add_argument("--low-score-cutoff", type=float, default=0.02, help="Raw-score norm cutoff for landmark filtering.")
    parser.add_argument("--bandwidth", type=int, default=100, help="Fourier bandwidth for each loop fit.")
    parser.add_argument("--regularization", type=float, default=5.0, help="Ridge regularization strength.")
    parser.add_argument("--curve-samples", type=int, default=1000, help="Samples per fitted curve in the overlay.")
    parser.add_argument("--noise", type=float, default=0.1, help="Gaussian noise added to the point cloud.")
    parser.add_argument("--persistence-landmarks", type=int, default=1000, help="n_perm passed to ripser for the persistence panel.")
    parser.add_argument("--skip-persistence", action="store_true", help="Skip the point-cloud/persistence SVG.")
    parser.add_argument("--skip-cup-product", action="store_true", help="Skip the cup-product table SVG.")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
