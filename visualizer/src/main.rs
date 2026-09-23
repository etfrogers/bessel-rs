use amos_bessel_rs::Scaling;
use amos_bessel_rs::amos::{
    complex_bessel_i_into, complex_bessel_j_into, complex_bessel_k_into, complex_bessel_y_into,
};
use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotBounds, PlotImage, PlotPoint, PlotPoints};
use num_complex::Complex64;
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, Ordering};

pub const CURVE_PALETTE: [egui::Color32; 8] = [
    egui::Color32::from_rgb(31, 119, 180),  // Tableau Blue
    egui::Color32::from_rgb(255, 127, 14),  // Tableau Orange
    egui::Color32::from_rgb(44, 160, 44),   // Tableau Green
    egui::Color32::from_rgb(214, 39, 40),   // Tableau Red
    egui::Color32::from_rgb(148, 103, 189), // Tableau Purple
    egui::Color32::from_rgb(140, 86, 75),   // Tableau Brown
    egui::Color32::from_rgb(227, 119, 194), // Tableau Pink
    egui::Color32::from_rgb(23, 190, 207),  // Tableau Cyan
];

pub fn get_curve_color(birth_id: i64) -> egui::Color32 {
    let n = CURVE_PALETTE.len() as i64;
    let idx = birth_id.rem_euclid(n) as usize;
    CURVE_PALETTE[idx]
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum ViewTab {
    Heatmap2D,
    Curves1D,
    Split,
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum BesselKind {
    J, // First kind
    Y, // Second kind
    I, // Modified first kind
    K, // Modified second kind
}

impl BesselKind {
    pub fn name(&self) -> &'static str {
        match self {
            Self::J => "J",
            Self::Y => "Y",
            Self::I => "I",
            Self::K => "K",
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            Self::J => "J_v (First Kind)",
            Self::Y => "Y_v (Second Kind)",
            Self::I => "I_v (Mod. 1st Kind)",
            Self::K => "K_v (Mod. 2nd Kind)",
        }
    }

    pub fn scaling_factor_2d(&self) -> &'static str {
        match self {
            Self::J => "e^{-|Im(z)|} · J_v(z)",
            Self::Y => "e^{-|Im(z)|} · Y_v(z)",
            Self::I => "e^{-|Re(z)|} · I_v(z)",
            Self::K => "e^z · K_v(z)",
        }
    }

    pub fn scaling_factor_1d(&self) -> &'static str {
        match self {
            Self::J => "e^{-|Im(x)|} · J_v(x) (= J_v(x))",
            Self::Y => "e^{-|Im(x)|} · Y_v(x) (= Y_v(x))",
            Self::I => "e^{-|x|} · I_v(x)",
            Self::K => "e^x · K_v(x)",
        }
    }
}

pub fn eval_bessel_sequence(
    kind: BesselKind,
    order: f64,
    z: Complex64,
    scaling: Scaling,
    out: &mut [Complex64],
) -> Result<(), amos_bessel_rs::BesselError<f64>> {
    match kind {
        BesselKind::J => complex_bessel_j_into(z, order, scaling, out)?,
        BesselKind::Y => complex_bessel_y_into(z, order, scaling, out)?,
        BesselKind::I => complex_bessel_i_into(z, order, scaling, out)?,
        BesselKind::K => complex_bessel_k_into(z, order, scaling, out)?,
    };
    Ok(())
}

pub fn eval_bessel_single(
    kind: BesselKind,
    order: f64,
    z: Complex64,
    scaling: Scaling,
) -> Result<Complex64, amos_bessel_rs::BesselError<f64>> {
    let mut buf = [Complex64::new(0.0, 0.0); 1];
    eval_bessel_sequence(kind, order, z, scaling, &mut buf)?;
    Ok(buf[0])
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum DomainMode {
    ComplexPlane, // x = Re(z), y = Im(z), fixed order
    RealVsOrder,  // x = Real(z), y = order v
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum Component {
    Magnitude,
    LogMagnitude,
    Phase,
    Real,
    Imag,
}

impl Component {
    pub fn extract(&self, val: Complex64) -> f64 {
        let v = match self {
            Component::Magnitude => val.norm(),
            Component::LogMagnitude => (val.norm() + 1e-12).ln(),
            Component::Phase => val.arg(),
            Component::Real => val.re,
            Component::Imag => val.im,
        };
        if v.is_finite() { v } else { f64::NAN }
    }
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum ColormapPreset {
    Turbo,
    Viridis,
    Plasma,
    Magma,
    Rainbow,
}

pub struct RawGrid {
    pub width: usize,
    pub height: usize,
    pub data: Vec<f64>,
}

pub struct BesselViewer {
    pub view_tab: ViewTab,

    // Animation state
    pub animating: bool,
    pub anim_speed: f64,          // cycles per second
    pub anim_2d_range: [f64; 2],  // [min_order, max_order] for 2D animation
    pub anim_2d_ping_pong: bool,  // ping-pong bounce vs forward loop
    pub anim_2d_forward: bool,    // direction indicator for ping-pong
    pub curve_anim_cycle: i64,    // completed integer cycles for rolling color tracking
    pub curve_anim_fraction: f64, // fractional order offset [0.0, 1.0) for 1D curves

    // 2D parameters
    pub kind: BesselKind,
    pub scaling: Scaling,
    pub domain_mode: DomainMode,
    pub component: Component,
    pub colormap: ColormapPreset,
    pub order: f64,
    pub x_range: [f64; 2],
    pub y_range: [f64; 2],
    pub val_range: [f64; 2],
    pub auto_range: bool,
    pub resolution: usize,
    pub dirty: bool,
    pub texture: Option<egui::TextureHandle>,
    pub gradient: Box<dyn colorgrad::Gradient + Send + Sync>,

    // 1D curve parameters
    pub curve_kind: BesselKind,
    pub curve_scaling: Scaling,
    pub curve_x_range: [f64; 2],
    pub curve_y_range: [f64; 2],
    pub lock_curve_y: bool,
    pub curve_reset_view: AtomicBool,
    pub curve_samples: usize,
    pub enabled_orders: Vec<(f64, bool)>, // list of (order, is_enabled)
    pub custom_order: f64,
    pub custom_order_enabled: bool,
}

pub fn build_gradient(preset: ColormapPreset) -> Box<dyn colorgrad::Gradient + Send + Sync> {
    match preset {
        ColormapPreset::Turbo => Box::new(colorgrad::preset::turbo()),
        ColormapPreset::Viridis => Box::new(colorgrad::preset::viridis()),
        ColormapPreset::Plasma => Box::new(colorgrad::preset::plasma()),
        ColormapPreset::Magma => Box::new(colorgrad::preset::magma()),
        ColormapPreset::Rainbow => Box::new(colorgrad::preset::rainbow()),
    }
}

impl BesselViewer {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let mut app = Self {
            view_tab: ViewTab::Heatmap2D,

            animating: false,
            anim_speed: 0.25, // 0.25 cycles/sec -> 4 seconds per cycle
            anim_2d_range: [0.0, 8.0],
            anim_2d_ping_pong: true,
            anim_2d_forward: true,
            curve_anim_cycle: 0,
            curve_anim_fraction: 0.0,

            kind: BesselKind::J,
            scaling: Scaling::Unscaled,
            domain_mode: DomainMode::ComplexPlane,
            component: Component::Magnitude,
            colormap: ColormapPreset::Turbo,
            order: 0.0,
            x_range: [-15.0, 15.0],
            y_range: [-15.0, 15.0],
            val_range: [0.0, 2.0],
            auto_range: false,
            resolution: 300,
            dirty: true,
            texture: None,
            gradient: build_gradient(ColormapPreset::Turbo),

            curve_kind: BesselKind::J,
            curve_scaling: Scaling::Unscaled,
            curve_x_range: [0.0, 20.0],
            curve_y_range: [-0.6, 1.1],
            lock_curve_y: true,
            curve_reset_view: AtomicBool::new(true),
            curve_samples: 600,
            enabled_orders: vec![
                (0.0, true),
                (1.0, true),
                (2.0, true),
                (3.0, true),
                (4.0, false),
                (5.0, false),
            ],
            custom_order: 0.5,
            custom_order_enabled: false,
        };
        app.update_texture(&cc.egui_ctx);
        app
    }

    pub fn step_animation(&mut self, dt: f64) {
        if !self.animating {
            return;
        }

        // 1. Advance 1D fractional order offset and cycle count for continuous color rolling
        let new_total =
            (self.curve_anim_cycle as f64 + self.curve_anim_fraction) + dt * self.anim_speed;
        self.curve_anim_cycle = new_total.floor() as i64;
        self.curve_anim_fraction = new_total.rem_euclid(1.0);

        // 2. Advance 2D complex plane order
        if (self.view_tab == ViewTab::Heatmap2D || self.view_tab == ViewTab::Split)
            && self.domain_mode == DomainMode::ComplexPlane
        {
            let low = self.anim_2d_range[0].min(self.anim_2d_range[1]);
            let high = self.anim_2d_range[0].max(self.anim_2d_range[1]);
            let span = high - low;

            if span > 1e-6 {
                let delta = dt * self.anim_speed * span;
                if self.anim_2d_ping_pong {
                    if self.anim_2d_forward {
                        self.order += delta;
                        if self.order >= high {
                            self.order = high;
                            self.anim_2d_forward = false;
                        }
                    } else {
                        self.order -= delta;
                        if self.order <= low {
                            self.order = low;
                            self.anim_2d_forward = true;
                        }
                    }
                } else {
                    self.order += delta;
                    if self.order > high {
                        self.order = low + (self.order - low).rem_euclid(span);
                    }
                }
                self.dirty = true;
            }
        }
    }

    pub fn eval_bessel(&self, x: f64, y: f64) -> f64 {
        let (order, z) = match self.domain_mode {
            DomainMode::ComplexPlane => (self.order, Complex64::new(x, y)),
            DomainMode::RealVsOrder => (y, Complex64::new(x, 0.0)),
        };

        let result = eval_bessel_single(self.kind, order, z, self.scaling);

        match result {
            Ok(val) => self.component.extract(val),
            Err(_) => f64::NAN,
        }
    }

    pub fn eval_1d_bessel(kind: BesselKind, order: f64, x: f64, scaling: Scaling) -> Option<f64> {
        let z = Complex64::new(x, 0.0);
        let res = eval_bessel_single(kind, order, z, scaling);

        match res {
            Ok(val) => {
                if val.im.abs() <= 1e-10 * (val.re.abs() + 1.0) && val.re.is_finite() {
                    Some(val.re)
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    pub fn compute_raw_grid(&self) -> RawGrid {
        let (x_min, x_max) = (self.x_range[0], self.x_range[1]);
        let (y_min, y_max) = (self.y_range[0], self.y_range[1]);
        let w = self.resolution;

        match self.domain_mode {
            DomainMode::ComplexPlane => {
                let h = self.resolution;
                let data = (0..h)
                    .into_par_iter()
                    .flat_map(|row| {
                        let v = y_max - (row as f64 / (h - 1) as f64) * (y_max - y_min);
                        (0..w)
                            .map(|col| {
                                let u = x_min + (col as f64 / (w - 1) as f64) * (x_max - x_min);
                                self.eval_bessel(u, v)
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect();
                RawGrid {
                    width: w,
                    height: h,
                    data,
                }
            }
            DomainMode::RealVsOrder => {
                let y_start = y_min.min(y_max);
                let y_end = y_min.max(y_max);
                let span = y_end - y_start;

                if span < 1e-6 || w < 2 {
                    let h = self.resolution;
                    let data = (0..h)
                        .into_par_iter()
                        .flat_map(|row| {
                            let v = y_max - (row as f64 / (h - 1) as f64) * (y_max - y_min);
                            (0..w)
                                .map(|col| {
                                    let u = x_min + (col as f64 / (w - 1) as f64) * (x_max - x_min);
                                    self.eval_bessel(u, v)
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect();
                    return RawGrid {
                        width: w,
                        height: h,
                        data,
                    };
                }

                // Discretize the unit interval into M fractional steps per unit order
                let m_sub = ((w - 1) as f64 / span).round().max(1.0) as usize;
                let h = (span * m_sub as f64).round() as usize + 1;

                // Compute columns in parallel across real z
                let cols: Vec<Vec<f64>> = (0..w)
                    .into_par_iter()
                    .map(|col| {
                        let u = x_min + (col as f64 / (w - 1) as f64) * (x_max - x_min);
                        let z = Complex64::new(u, 0.0);
                        let mut col_vals = vec![f64::NAN; h];

                        for m in 0..m_sub {
                            let k_count = (h - 1 - m) / m_sub + 1;
                            let nu_start = y_start + (m as f64 / m_sub as f64);
                            let mut buf = vec![Complex64::new(0.0, 0.0); k_count];

                            if eval_bessel_sequence(self.kind, nu_start, z, self.scaling, &mut buf)
                                .is_ok()
                            {
                                for (k, &c_val) in buf.iter().enumerate() {
                                    let r = m + k * m_sub;
                                    if r < h {
                                        col_vals[r] = self.component.extract(c_val);
                                    }
                                }
                            }
                        }
                        col_vals
                    })
                    .collect();

                // Reconstruct 2D row-major raster
                // Texture row 0 corresponds to y_max, row h-1 corresponds to y_min
                let mut data = vec![f64::NAN; w * h];
                for img_row in 0..h {
                    let r = if y_max >= y_min {
                        (h - 1) - img_row
                    } else {
                        img_row
                    };
                    for col in 0..w {
                        data[img_row * w + col] = cols[col][r];
                    }
                }

                RawGrid {
                    width: w,
                    height: h,
                    data,
                }
            }
        }
    }

    fn update_texture(&mut self, ctx: &egui::Context) {
        let grid = self.compute_raw_grid();
        let (width, height) = (grid.width, grid.height);
        let raw_values = grid.data;

        // Normalization & clipping bounds
        let (val_min, val_max) = if self.auto_range {
            let mut valid: Vec<f64> = raw_values
                .iter()
                .copied()
                .filter(|v| v.is_finite())
                .collect();
            if valid.is_empty() {
                (0.0, 1.0)
            } else {
                valid.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let p2 = valid[(valid.len() as f64 * 0.02) as usize];
                let p98 = valid[(valid.len() as f64 * 0.98).min((valid.len() - 1) as f64) as usize];
                if (p98 - p2).abs() < 1e-6 {
                    (p2 - 0.5, p2 + 0.5)
                } else {
                    (p2, p98)
                }
            }
        } else {
            (self.val_range[0], self.val_range[1])
        };

        let range_span = if (val_max - val_min).abs() < 1e-9 {
            1.0
        } else {
            val_max - val_min
        };

        let mut pixels = Vec::with_capacity(width * height * 4);
        for &val in &raw_values {
            if val.is_nan() {
                pixels.extend_from_slice(&[30, 30, 30, 255]);
            } else {
                let t = ((val - val_min) / range_span).clamp(0.0, 1.0);
                let rgba = self.gradient.at(t as f32).to_rgba8();
                pixels.extend_from_slice(&rgba);
            }
        }

        let image = egui::ColorImage::from_rgba_unmultiplied([width, height], &pixels);
        self.texture = Some(ctx.load_texture("heatmap", image, egui::TextureOptions::LINEAR));
        self.dirty = false;
    }

    fn show_heatmap_plot(&self, ui: &mut egui::Ui) {
        let x_label = match self.domain_mode {
            DomainMode::ComplexPlane => "Re(z)",
            DomainMode::RealVsOrder => "Real(z)",
        };
        let y_label = match self.domain_mode {
            DomainMode::ComplexPlane => "Im(z)",
            DomainMode::RealVsOrder => "Order (v)",
        };

        let plot = Plot::new("bessel_surface")
            .data_aspect(1.0)
            .show_axes(true)
            .show_grid(true)
            .x_axis_label(x_label)
            .y_axis_label(y_label);

        plot.show(ui, |plot_ui| {
            if let Some(ref texture) = self.texture {
                let center_x = (self.x_range[0] + self.x_range[1]) * 0.5;
                let center_y = (self.y_range[0] + self.y_range[1]) * 0.5;
                let width = (self.x_range[1] - self.x_range[0]).abs();
                let height = (self.y_range[1] - self.y_range[0]).abs();

                plot_ui.image(PlotImage::new(
                    "bessel_heatmap",
                    texture.id(),
                    PlotPoint::new(center_x, center_y),
                    egui::vec2(width as f32, height as f32),
                ));
            }
        });
    }

    fn show_curves_plot(&self, ui: &mut egui::Ui) {
        let y_axis_label = if self.curve_scaling == Scaling::Scaled {
            match self.curve_kind {
                BesselKind::J => "Scaled J_v(x) [e^{-|Im(x)|} J_v(x)]",
                BesselKind::Y => "Scaled Y_v(x) [e^{-|Im(x)|} Y_v(x)]",
                BesselKind::I => "Scaled I_v(x) [e^{-|x|} I_v(x)]",
                BesselKind::K => "Scaled K_v(x) [e^x K_v(x)]",
            }
        } else {
            match self.curve_kind {
                BesselKind::J => "J_v(x)",
                BesselKind::Y => "Y_v(x)",
                BesselKind::I => "I_v(x)",
                BesselKind::K => "K_v(x)",
            }
        };

        let plot = Plot::new("bessel_curves")
            .legend(Legend::default())
            .show_axes(true)
            .show_grid(true)
            .x_axis_label("x")
            .y_axis_label(y_axis_label);

        plot.show(ui, |plot_ui| {
            if self.lock_curve_y {
                let cur = plot_ui.plot_bounds();
                let is_default_bounds =
                    (cur.min()[0] - (-1.0)).abs() < 1e-6 && (cur.max()[0] - 1.0).abs() < 1e-6;
                let reset =
                    self.curve_reset_view.swap(false, Ordering::Relaxed) || is_default_bounds;

                let min_x = if reset {
                    self.curve_x_range[0]
                } else {
                    cur.min()[0]
                };
                let max_x = if reset {
                    self.curve_x_range[1]
                } else {
                    cur.max()[0]
                };
                let min_y = self.curve_y_range[0].min(self.curve_y_range[1]);
                let max_y = self.curve_y_range[0].max(self.curve_y_range[1]);
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([min_x, min_y], [max_x, max_y]));
            }

            let n_pts = self.curve_samples.max(50);
            let (x_min, x_max) = (self.curve_x_range[0], self.curve_x_range[1]);

            struct ActiveOrder {
                order: f64,
                color: egui::Color32,
            }

            let mut active_orders: Vec<ActiveOrder> = Vec::new();
            for &(order, enabled) in &self.enabled_orders {
                if enabled {
                    // Continuous rolling color: birth_id tracks the wave's continuous propagation
                    let birth_id = self.curve_anim_cycle - order.round() as i64;
                    let color = get_curve_color(birth_id);
                    active_orders.push(ActiveOrder {
                        order: order + self.curve_anim_fraction,
                        color,
                    });
                }
            }
            if self.custom_order_enabled {
                let birth_id = self.curve_anim_cycle - self.custom_order.round() as i64;
                let color = get_curve_color(birth_id);
                active_orders.push(ActiveOrder {
                    order: self.custom_order + self.curve_anim_fraction,
                    color,
                });
            }

            active_orders.sort_by(|a, b| {
                a.order
                    .partial_cmp(&b.order)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            active_orders.dedup_by(|a, b| (a.order - b.order).abs() < 1e-9);

            // Group active orders into sequence batches separated by integer steps of 1.0
            struct SequenceBatch {
                base_order: f64,
                seq_len: usize,
                orders: Vec<(f64, usize, egui::Color32)>, // (order, buffer_offset, color)
            }

            let mut batches: Vec<SequenceBatch> = Vec::new();
            for item in &active_orders {
                let mut matched = false;
                for batch in &mut batches {
                    let diff = item.order - batch.base_order;
                    let rounded = diff.round();
                    if (diff - rounded).abs() < 1e-6 && rounded >= 0.0 {
                        let offset = rounded as usize;
                        batch.orders.push((item.order, offset, item.color));
                        if offset + 1 > batch.seq_len {
                            batch.seq_len = offset + 1;
                        }
                        matched = true;
                        break;
                    }
                }
                if !matched {
                    batches.push(SequenceBatch {
                        base_order: item.order,
                        seq_len: 1,
                        orders: vec![(item.order, 0, item.color)],
                    });
                }
            }

            // Evaluate each batch over all x samples in single sequence calls
            for batch in batches {
                let mut batch_points: Vec<Vec<[f64; 2]>> =
                    vec![Vec::with_capacity(n_pts); batch.orders.len()];
                let mut buf = vec![Complex64::new(0.0, 0.0); batch.seq_len];

                for i in 0..n_pts {
                    let t = i as f64 / (n_pts - 1) as f64;
                    let x = x_min + t * (x_max - x_min);
                    let z = Complex64::new(x, 0.0);

                    if eval_bessel_sequence(
                        self.curve_kind,
                        batch.base_order,
                        z,
                        self.curve_scaling,
                        &mut buf,
                    )
                    .is_ok()
                    {
                        for (ord_idx, &(_order, offset, _color)) in batch.orders.iter().enumerate()
                        {
                            let val = buf[offset];
                            if val.im.abs() <= 1e-10 * (val.re.abs() + 1.0) && val.re.is_finite() {
                                batch_points[ord_idx].push([x, val.re]);
                            }
                        }
                    }
                }

                for (ord_idx, &(order, _, color)) in batch.orders.iter().enumerate() {
                    let points = PlotPoints::new(std::mem::take(&mut batch_points[ord_idx]));
                    let label = if (order.fract()).abs() < 1e-6 {
                        if self.curve_scaling == Scaling::Scaled {
                            format!("{}_{{{:.0}}}(x) [scaled]", self.curve_kind.name(), order)
                        } else {
                            format!("{}_{{{:.0}}}(x)", self.curve_kind.name(), order)
                        }
                    } else if self.curve_scaling == Scaling::Scaled {
                        format!("{}_{{{:.2}}}(x) [scaled]", self.curve_kind.name(), order)
                    } else {
                        format!("{}_{{{:.2}}}(x)", self.curve_kind.name(), order)
                    };

                    plot_ui.line(Line::new(label, points).color(color));
                }
            }
        });
    }
}

impl eframe::App for BesselViewer {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        let dt = ui.input(|i| i.stable_dt).min(0.1) as f64;

        if self.animating {
            self.step_animation(dt);
            ui.ctx().request_repaint();
        }

        // Left Controls Panel
        egui::Panel::left("controls")
            .resizable(true)
            .default_size(320.0)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("Bessel Visualizer");
                    ui.add_space(6.0);

                    // View Mode Tabs
                    ui.horizontal(|ui| {
                        ui.selectable_value(&mut self.view_tab, ViewTab::Heatmap2D, "2D Heatmap");
                        ui.selectable_value(&mut self.view_tab, ViewTab::Curves1D, "1D Curves");
                        ui.selectable_value(&mut self.view_tab, ViewTab::Split, "Split View");
                    });

                    ui.add_space(8.0);
                    ui.separator();

                    // Section 0: Animation Controls
                    ui.heading("Animation Controls");
                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        let play_pause_label = if self.animating {
                            "⏸ Pause"
                        } else {
                            "▶ Play"
                        };
                        if ui
                            .button(egui::RichText::new(play_pause_label).strong())
                            .clicked()
                        {
                            self.animating = !self.animating;
                        }
                        if ui.button("↺ Reset").clicked() {
                            self.curve_anim_fraction = 0.0;
                            self.curve_anim_cycle = 0;
                            self.order = self.anim_2d_range[0];
                            self.anim_2d_forward = true;
                            self.curve_reset_view.store(true, Ordering::Relaxed);
                            self.dirty = true;
                        }
                    });

                    ui.add(
                        egui::Slider::new(&mut self.anim_speed, 0.02..=2.0).text("Speed (cyc/s)"),
                    );

                    if self.view_tab == ViewTab::Curves1D || self.view_tab == ViewTab::Split {
                        ui.horizontal(|ui| {
                            ui.label("1D Fraction (α):");
                            ui.add(
                                egui::Slider::new(&mut self.curve_anim_fraction, 0.0..=0.999)
                                    .text("")
                                    .show_value(true),
                            );
                        });
                        ui.checkbox(&mut self.lock_curve_y, "Lock 1D Y scale");
                    }

                    if (self.view_tab == ViewTab::Heatmap2D || self.view_tab == ViewTab::Split)
                        && self.domain_mode == DomainMode::ComplexPlane
                    {
                        ui.collapsing("2D Order Cycling", |ui| {
                            ui.horizontal(|ui| {
                                ui.label("Min:");
                                ui.add(egui::DragValue::new(&mut self.anim_2d_range[0]).speed(0.1));
                                ui.label("Max:");
                                ui.add(egui::DragValue::new(&mut self.anim_2d_range[1]).speed(0.1));
                            });
                            ui.checkbox(&mut self.anim_2d_ping_pong, "Ping-Pong (bounce)");
                        });
                    }

                    ui.add_space(8.0);
                    ui.separator();

                    // Section 1: 1D Curve Settings
                    if self.view_tab == ViewTab::Curves1D || self.view_tab == ViewTab::Split {
                        ui.heading("1D Curves Settings");
                        ui.add_space(4.0);

                        egui::ComboBox::from_label("Function Kind")
                            .selected_text(self.curve_kind.label())
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut self.curve_kind,
                                    BesselKind::J,
                                    BesselKind::J.label(),
                                );
                                ui.selectable_value(
                                    &mut self.curve_kind,
                                    BesselKind::Y,
                                    BesselKind::Y.label(),
                                );
                                ui.selectable_value(
                                    &mut self.curve_kind,
                                    BesselKind::I,
                                    BesselKind::I.label(),
                                );
                                ui.selectable_value(
                                    &mut self.curve_kind,
                                    BesselKind::K,
                                    BesselKind::K.label(),
                                );
                            });

                        let mut curve_scaled = self.curve_scaling == Scaling::Scaled;
                        if ui
                            .checkbox(&mut curve_scaled, "Exponential scaling")
                            .on_hover_text(format!(
                                "Scale to remove exponential growth/decay:\n{}",
                                self.curve_kind.scaling_factor_1d()
                            ))
                            .changed()
                        {
                            self.curve_scaling = if curve_scaled {
                                Scaling::Scaled
                            } else {
                                Scaling::Unscaled
                            };
                        }
                        if curve_scaled {
                            ui.label(
                                egui::RichText::new(format!(
                                    "Scaling factor: {}",
                                    self.curve_kind.scaling_factor_1d()
                                ))
                                .weak()
                                .small(),
                            );
                        }

                        ui.add_space(4.0);
                        ui.label("Orders to Overlay (Batch Evaluated):");
                        ui.horizontal_wrapped(|ui| {
                            for (order, enabled) in &mut self.enabled_orders {
                                let birth_id = self.curve_anim_cycle - order.round() as i64;
                                let color = get_curve_color(birth_id);
                                ui.horizontal(|ui| {
                                    ui.colored_label(color, "■");
                                    ui.checkbox(enabled, format!("v = {:.0}", order));
                                });
                            }
                        });

                        ui.horizontal(|ui| {
                            let birth_id = self.curve_anim_cycle - self.custom_order.round() as i64;
                            let custom_color = get_curve_color(birth_id);
                            ui.colored_label(custom_color, "■");
                            ui.checkbox(&mut self.custom_order_enabled, "Custom v:");
                            if self.custom_order_enabled {
                                ui.add(egui::DragValue::new(&mut self.custom_order).speed(0.05));
                            }
                        });

                        ui.add_space(4.0);
                        if ui
                            .add(
                                egui::Slider::new(&mut self.curve_x_range[0], -50.0..=10.0)
                                    .text("Curve X Min"),
                            )
                            .changed()
                        {
                            self.curve_reset_view.store(true, Ordering::Relaxed);
                        }
                        if ui
                            .add(
                                egui::Slider::new(&mut self.curve_x_range[1], 1.0..=100.0)
                                    .text("Curve X Max"),
                            )
                            .changed()
                        {
                            self.curve_reset_view.store(true, Ordering::Relaxed);
                        }
                        ui.add(
                            egui::Slider::new(&mut self.curve_samples, 100..=2000).text("Samples"),
                        );

                        ui.add_space(4.0);
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut self.lock_curve_y, "Lock Y scale");
                            if ui.small_button("Fit View").clicked() {
                                self.curve_reset_view.store(true, Ordering::Relaxed);
                            }
                        });
                        if self.lock_curve_y {
                            let cur_y_min = self.curve_y_range[0];
                            let cur_y_max = self.curve_y_range[1];
                            ui.horizontal(|ui| {
                                ui.label("Y Min:");
                                ui.add(
                                    egui::DragValue::new(&mut self.curve_y_range[0])
                                        .speed(0.05)
                                        .range(-1000.0..=cur_y_max - 0.01),
                                );
                                ui.label("Y Max:");
                                ui.add(
                                    egui::DragValue::new(&mut self.curve_y_range[1])
                                        .speed(0.05)
                                        .range(cur_y_min + 0.01..=1000.0),
                                );
                            });

                            ui.horizontal(|ui| {
                                ui.label("Presets:");
                                if ui.small_button("[-0.6, 1.1]").clicked() {
                                    self.curve_y_range = [-0.6, 1.1];
                                }
                                if ui.small_button("[-1.5, 1.5]").clicked() {
                                    self.curve_y_range = [-1.5, 1.5];
                                }
                                if ui.small_button("[-5.0, 5.0]").clicked() {
                                    self.curve_y_range = [-5.0, 5.0];
                                }
                            });
                        }

                        ui.add_space(8.0);
                        ui.separator();
                    }

                    // Section 2: 2D Heatmap Settings
                    if self.view_tab == ViewTab::Heatmap2D || self.view_tab == ViewTab::Split {
                        ui.heading("2D Heatmap Settings");
                        ui.add_space(4.0);

                        egui::ComboBox::from_label("2D Type")
                            .selected_text(self.kind.label())
                            .show_ui(ui, |ui| {
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.kind,
                                        BesselKind::J,
                                        BesselKind::J.label(),
                                    )
                                    .changed();
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.kind,
                                        BesselKind::Y,
                                        BesselKind::Y.label(),
                                    )
                                    .changed();
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.kind,
                                        BesselKind::I,
                                        BesselKind::I.label(),
                                    )
                                    .changed();
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.kind,
                                        BesselKind::K,
                                        BesselKind::K.label(),
                                    )
                                    .changed();
                            });

                        let mut scaled_2d = self.scaling == Scaling::Scaled;
                        if ui
                            .checkbox(&mut scaled_2d, "Exponential scaling")
                            .on_hover_text(format!(
                                "Scale to remove exponential growth/decay:\n{}",
                                self.kind.scaling_factor_2d()
                            ))
                            .changed()
                        {
                            self.scaling = if scaled_2d {
                                Scaling::Scaled
                            } else {
                                Scaling::Unscaled
                            };
                            self.dirty = true;
                        }
                        if scaled_2d {
                            ui.label(
                                egui::RichText::new(format!(
                                    "Scaling factor: {}",
                                    self.kind.scaling_factor_2d()
                                ))
                                .weak()
                                .small(),
                            );
                        }

                        ui.add_space(4.0);
                        egui::ComboBox::from_label("2D Domain")
                            .selected_text(match self.domain_mode {
                                DomainMode::ComplexPlane => "Complex Plane: z = x + iy",
                                DomainMode::RealVsOrder => "Real vs Order: z = x, v = y",
                            })
                            .show_ui(ui, |ui| {
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.domain_mode,
                                        DomainMode::ComplexPlane,
                                        "Complex Plane: z = x + iy",
                                    )
                                    .changed();
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.domain_mode,
                                        DomainMode::RealVsOrder,
                                        "Real vs Order: z = x, v = y",
                                    )
                                    .changed();
                            });

                        if self.domain_mode == DomainMode::ComplexPlane {
                            self.dirty |= ui
                                .add(
                                    egui::Slider::new(&mut self.order, -10.0..=10.0)
                                        .text("Order (v)"),
                                )
                                .changed();
                        }

                        egui::ComboBox::from_label("Component")
                            .selected_text(match self.component {
                                Component::Magnitude => "Magnitude |f(z)|",
                                Component::LogMagnitude => "Log Magnitude ln|f(z)|",
                                Component::Phase => "Phase arg(f(z))",
                                Component::Real => "Real Part Re(f(z))",
                                Component::Imag => "Imag Part Im(f(z))",
                            })
                            .show_ui(ui, |ui| {
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.component,
                                        Component::Magnitude,
                                        "Magnitude |f(z)|",
                                    )
                                    .changed();
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.component,
                                        Component::LogMagnitude,
                                        "Log Magnitude ln|f(z)|",
                                    )
                                    .changed();
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.component,
                                        Component::Phase,
                                        "Phase arg(f(z))",
                                    )
                                    .changed();
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.component,
                                        Component::Real,
                                        "Real Part Re(f(z))",
                                    )
                                    .changed();
                                self.dirty |= ui
                                    .selectable_value(
                                        &mut self.component,
                                        Component::Imag,
                                        "Imag Part Im(f(z))",
                                    )
                                    .changed();
                            });

                        ui.add_space(6.0);

                        let mut colormap_changed = false;
                        egui::ComboBox::from_label("Colormap")
                            .selected_text(format!("{:?}", self.colormap))
                            .show_ui(ui, |ui| {
                                colormap_changed |= ui
                                    .selectable_value(
                                        &mut self.colormap,
                                        ColormapPreset::Turbo,
                                        "Turbo",
                                    )
                                    .changed();
                                colormap_changed |= ui
                                    .selectable_value(
                                        &mut self.colormap,
                                        ColormapPreset::Viridis,
                                        "Viridis",
                                    )
                                    .changed();
                                colormap_changed |= ui
                                    .selectable_value(
                                        &mut self.colormap,
                                        ColormapPreset::Plasma,
                                        "Plasma",
                                    )
                                    .changed();
                                colormap_changed |= ui
                                    .selectable_value(
                                        &mut self.colormap,
                                        ColormapPreset::Magma,
                                        "Magma",
                                    )
                                    .changed();
                                colormap_changed |= ui
                                    .selectable_value(
                                        &mut self.colormap,
                                        ColormapPreset::Rainbow,
                                        "Rainbow",
                                    )
                                    .changed();
                            });

                        if colormap_changed {
                            self.gradient = build_gradient(self.colormap);
                            self.dirty = true;
                        }

                        self.dirty |= ui
                            .add(
                                egui::Slider::new(&mut self.resolution, 100..=800).text("Grid Res"),
                            )
                            .changed();

                        self.dirty |= ui
                            .checkbox(&mut self.auto_range, "Auto dynamic range (2%-98%)")
                            .changed();

                        if !self.auto_range {
                            self.dirty |= ui
                                .add(
                                    egui::Slider::new(&mut self.val_range[0], -10.0..=10.0)
                                        .text("Val Min"),
                                )
                                .changed();
                            self.dirty |= ui
                                .add(
                                    egui::Slider::new(&mut self.val_range[1], -10.0..=10.0)
                                        .text("Val Max"),
                                )
                                .changed();
                        }

                        ui.add_space(6.0);
                        ui.label("Domain Limits:");
                        self.dirty |= ui
                            .add(egui::Slider::new(&mut self.x_range[0], -50.0..=0.0).text("X Min"))
                            .changed();
                        self.dirty |= ui
                            .add(egui::Slider::new(&mut self.x_range[1], 0.0..=50.0).text("X Max"))
                            .changed();
                        self.dirty |= ui
                            .add(egui::Slider::new(&mut self.y_range[0], -50.0..=0.0).text("Y Min"))
                            .changed();
                        self.dirty |= ui
                            .add(egui::Slider::new(&mut self.y_range[1], 0.0..=50.0).text("Y Max"))
                            .changed();

                        if self.dirty {
                            self.update_texture(&ctx);
                        }
                    }
                });
            });

        // Central Area
        egui::CentralPanel::default().show(ui, |ui| match self.view_tab {
            ViewTab::Heatmap2D => {
                let scaled_str = if self.scaling == Scaling::Scaled {
                    " [Exponentially Scaled]"
                } else {
                    ""
                };
                ui.heading(format!(
                    "2D Heatmap Surface: {} (v = {:.2}){}",
                    self.kind.label(),
                    self.order,
                    scaled_str
                ));
                self.show_heatmap_plot(ui);
            }
            ViewTab::Curves1D => {
                let scaled_str = if self.curve_scaling == Scaling::Scaled {
                    " [Exponentially Scaled]"
                } else {
                    ""
                };
                let frac_str = if self.curve_anim_fraction > 1e-4 {
                    format!(" (+{:.2})", self.curve_anim_fraction)
                } else {
                    String::new()
                };
                ui.heading(format!(
                    "1D Overlaid Curves: {}{}{}",
                    self.curve_kind.label(),
                    frac_str,
                    scaled_str
                ));
                self.show_curves_plot(ui);
            }
            ViewTab::Split => {
                ui.columns(2, |columns| {
                    columns[0].vertical(|ui| {
                        let scaled_str = if self.scaling == Scaling::Scaled {
                            " [Scaled]"
                        } else {
                            ""
                        };
                        ui.heading(format!(
                            "2D Heatmap: {} (v = {:.2}){}",
                            self.kind.label(),
                            self.order,
                            scaled_str
                        ));
                        self.show_heatmap_plot(ui);
                    });
                    columns[1].vertical(|ui| {
                        let scaled_str = if self.curve_scaling == Scaling::Scaled {
                            " [Scaled]"
                        } else {
                            ""
                        };
                        let frac_str = if self.curve_anim_fraction > 1e-4 {
                            format!(" (+{:.2})", self.curve_anim_fraction)
                        } else {
                            String::new()
                        };
                        ui.heading(format!(
                            "1D Curves: {}{}{}",
                            self.curve_kind.label(),
                            frac_str,
                            scaled_str
                        ));
                        self.show_curves_plot(ui);
                    });
                });
            }
        });
    }
}

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_title("Bessel Function Visualizer"),
        ..Default::default()
    };
    eframe::run_native(
        "Bessel Function Visualizer",
        options,
        Box::new(|cc| Ok(Box::new(BesselViewer::new(cc)))),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_viewer() -> BesselViewer {
        BesselViewer {
            view_tab: ViewTab::Heatmap2D,
            animating: false,
            anim_speed: 1.0,
            anim_2d_range: [0.0, 5.0],
            anim_2d_ping_pong: true,
            anim_2d_forward: true,
            curve_anim_cycle: 0,
            curve_anim_fraction: 0.0,
            kind: BesselKind::J,
            scaling: Scaling::Unscaled,
            domain_mode: DomainMode::ComplexPlane,
            component: Component::Magnitude,
            colormap: ColormapPreset::Turbo,
            order: 0.0,
            x_range: [-5.0, 5.0],
            y_range: [-5.0, 5.0],
            val_range: [0.0, 1.0],
            auto_range: false,
            resolution: 10,
            dirty: false,
            texture: None,
            gradient: build_gradient(ColormapPreset::Turbo),
            curve_kind: BesselKind::J,
            curve_scaling: Scaling::Unscaled,
            curve_x_range: [0.0, 20.0],
            curve_y_range: [-0.6, 1.1],
            lock_curve_y: true,
            curve_reset_view: AtomicBool::new(true),
            curve_samples: 100,
            enabled_orders: vec![(0.0, true)],
            custom_order: 0.5,
            custom_order_enabled: false,
        }
    }

    #[test]
    fn test_bessel_j_origin() {
        let viewer = make_test_viewer();
        // J_0(0) = 1.0
        let val = viewer.eval_bessel(0.0, 0.0);
        assert!((val - 1.0).abs() < 1e-6, "J_0(0) should be 1.0, got {val}");
    }

    #[test]
    fn test_lock_curve_y_scale_default() {
        let viewer = make_test_viewer();
        assert!(viewer.lock_curve_y);
        assert_eq!(viewer.curve_y_range, [-0.6, 1.1]);
        assert_eq!(viewer.curve_x_range, [0.0, 20.0]);
        assert_eq!(viewer.curve_anim_cycle, 0);
        assert!(viewer.curve_reset_view.load(Ordering::Relaxed));
    }

    #[test]
    fn test_curve_palette_colors_distinct() {
        assert_eq!(CURVE_PALETTE.len(), 8);
        for i in 0..CURVE_PALETTE.len() {
            for j in (i + 1)..CURVE_PALETTE.len() {
                assert_ne!(
                    CURVE_PALETTE[i], CURVE_PALETTE[j],
                    "Palette colors {i} and {j} must be distinct"
                );
            }
        }
    }

    #[test]
    fn test_curve_color_rolling_continuity() {
        // Curve born at order 0 in cycle 0:
        let color_cycle0_order0 = get_curve_color(0 - 0); // birth_id = 0

        // In cycle 1, the curve that moved to order 1 should keep the exact same color!
        let color_cycle1_order1 = get_curve_color(1 - 1); // birth_id = 0
        assert_eq!(
            color_cycle0_order0, color_cycle1_order1,
            "The curve that was at order 0 must retain its color as it rolls on to order 1"
        );

        // In cycle 1, a brand new line is introduced at order 0 with a fresh color!
        let color_cycle1_order0 = get_curve_color(1 - 0); // birth_id = 1
        assert_ne!(
            color_cycle1_order0, color_cycle1_order1,
            "The new line entering at order 0 must have a distinct new color"
        );

        // In cycle 2, the original curve is now at order 2, still with its original color!
        let color_cycle2_order2 = get_curve_color(2 - 2); // birth_id = 0
        assert_eq!(color_cycle0_order0, color_cycle2_order2);
    }

    #[test]
    fn test_compute_raw_grid_complex_plane() {
        let mut viewer = make_test_viewer();
        viewer.resolution = 20;
        viewer.order = 1.0;
        let grid = viewer.compute_raw_grid();
        assert_eq!(grid.data.len(), 400);
        assert_eq!(grid.width, 20);
        assert_eq!(grid.height, 20);
        assert!(grid.data.iter().any(|&v| v.is_finite()));
    }

    #[test]
    fn test_compute_raw_grid_real_vs_order_sequence() {
        let mut viewer = make_test_viewer();
        viewer.domain_mode = DomainMode::RealVsOrder;
        viewer.component = Component::Real;
        viewer.x_range = [1.0, 10.0];
        viewer.y_range = [0.0, 5.0];
        viewer.resolution = 25;

        let grid = viewer.compute_raw_grid();
        assert_eq!(grid.data.len(), grid.width * grid.height);
        assert_eq!(grid.width, 25);
        assert!(grid.height > 0);
        assert!(grid.data.iter().all(|&v| v.is_finite()));

        // Sample check: J_0(2.0)
        let j0_2 =
            BesselViewer::eval_1d_bessel(BesselKind::J, 0.0, 2.0, Scaling::Unscaled).unwrap();
        let eval_direct = viewer.eval_bessel(2.0, 0.0);
        assert!((j0_2 - eval_direct).abs() < 1e-10);
    }

    #[test]
    fn test_compute_raw_grid_scaled() {
        let mut viewer = make_test_viewer();
        viewer.kind = BesselKind::I;
        viewer.scaling = Scaling::Scaled;
        viewer.resolution = 20;

        let grid = viewer.compute_raw_grid();
        assert_eq!(grid.data.len(), 400);
        assert!(grid.data.iter().all(|&v| v.is_finite() && !v.is_nan()));
    }

    #[test]
    fn test_eval_1d_bessel() {
        let j0_0 = BesselViewer::eval_1d_bessel(BesselKind::J, 0.0, 0.0, Scaling::Unscaled);
        assert_eq!(j0_0, Some(1.0));

        let j1_0 = BesselViewer::eval_1d_bessel(BesselKind::J, 1.0, 0.0, Scaling::Unscaled);
        assert!((j1_0.unwrap() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_eval_scaled_bessel_i_k() {
        let unscaled_i =
            BesselViewer::eval_1d_bessel(BesselKind::I, 0.0, 2.0, Scaling::Unscaled).unwrap();
        let scaled_i =
            BesselViewer::eval_1d_bessel(BesselKind::I, 0.0, 2.0, Scaling::Scaled).unwrap();
        let expected_scaled_i = unscaled_i * (-2.0_f64).exp();
        assert!(
            (scaled_i - expected_scaled_i).abs() < 1e-12,
            "Scaled I_0(2.0) mismatch: {scaled_i} vs expected {expected_scaled_i}"
        );

        let unscaled_k =
            BesselViewer::eval_1d_bessel(BesselKind::K, 0.0, 2.0, Scaling::Unscaled).unwrap();
        let scaled_k =
            BesselViewer::eval_1d_bessel(BesselKind::K, 0.0, 2.0, Scaling::Scaled).unwrap();
        let expected_scaled_k = unscaled_k * 2.0_f64.exp();
        assert!(
            (scaled_k - expected_scaled_k).abs() < 1e-12,
            "Scaled K_0(2.0) mismatch: {scaled_k} vs expected {expected_scaled_k}"
        );
    }

    #[test]
    fn test_sequence_batching_equivalence() {
        let z = Complex64::new(3.5, 0.0);
        let mut seq_buf = [Complex64::new(0.0, 0.0); 4];
        eval_bessel_sequence(BesselKind::J, 0.0, z, Scaling::Unscaled, &mut seq_buf).unwrap();

        for (k, expected_val) in seq_buf.iter().enumerate() {
            let single_val =
                eval_bessel_single(BesselKind::J, k as f64, z, Scaling::Unscaled).unwrap();
            assert!(
                (expected_val - single_val).norm() < 1e-12,
                "Order {k} sequence mismatch: {expected_val} vs {single_val}"
            );
        }
    }

    #[test]
    fn test_animation_step_1d_and_2d() {
        let mut viewer = make_test_viewer();
        viewer.animating = true;
        viewer.anim_speed = 1.0;
        viewer.anim_2d_range = [0.0, 4.0];
        viewer.anim_2d_ping_pong = true;
        viewer.anim_2d_forward = true;
        viewer.curve_anim_cycle = 0;
        viewer.curve_anim_fraction = 0.0;
        viewer.order = 0.0;

        // Step by 0.25 seconds:
        // 1D fractional offset should advance to 0.25, cycle stays 0
        // 2D order should advance: 0.0 + 0.25 * 1.0 * 4.0 = 1.0
        viewer.step_animation(0.25);
        assert!((viewer.curve_anim_fraction - 0.25).abs() < 1e-9);
        assert_eq!(viewer.curve_anim_cycle, 0);
        assert!((viewer.order - 1.0).abs() < 1e-9);
        assert!(viewer.dirty);

        // Step past 1.0 seconds total (dt = 0.8):
        // Total progress: 0.25 + 0.8 = 1.05 -> cycle 1, fraction 0.05
        viewer.step_animation(0.8);
        assert_eq!(viewer.curve_anim_cycle, 1);
        assert!((viewer.curve_anim_fraction - 0.05).abs() < 1e-9);

        // Test 2D ping pong boundary reversal:
        // Current order: 1.0 + 0.8 * 4.0 = 4.2 -> clamps to 4.0 and reverses forward flag
        assert!((viewer.order - 4.0).abs() < 1e-9);
        assert!(!viewer.anim_2d_forward);

        // Next step should decrease order:
        viewer.step_animation(0.25);
        assert!((viewer.order - 3.0).abs() < 1e-9);
    }
}
