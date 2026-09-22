use amos_bessel_rs::Scaling;
use amos_bessel_rs::amos::{
    complex_bessel_i_into, complex_bessel_j_into, complex_bessel_k_into, complex_bessel_y_into,
};
use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotImage, PlotPoint, PlotPoints};
use num_complex::Complex64;
use rayon::prelude::*;

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

pub fn eval_bessel_single(
    kind: BesselKind,
    order: f64,
    z: Complex64,
    scaling: Scaling,
) -> Result<Complex64, amos_bessel_rs::BesselError<f64>> {
    let mut buf = [Complex64::new(0.0, 0.0); 1];
    match kind {
        BesselKind::J => complex_bessel_j_into(z, order, scaling, &mut buf)?,
        BesselKind::Y => complex_bessel_y_into(z, order, scaling, &mut buf)?,
        BesselKind::I => complex_bessel_i_into(z, order, scaling, &mut buf)?,
        BesselKind::K => complex_bessel_k_into(z, order, scaling, &mut buf)?,
    };
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

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum ColormapPreset {
    Turbo,
    Viridis,
    Plasma,
    Magma,
    Rainbow,
}

pub struct BesselViewer {
    view_tab: ViewTab,
    // 2D parameters
    kind: BesselKind,
    scaling: Scaling,
    domain_mode: DomainMode,
    component: Component,
    colormap: ColormapPreset,
    order: f64,
    x_range: [f64; 2],
    y_range: [f64; 2],
    val_range: [f64; 2],
    auto_range: bool,
    resolution: usize,
    dirty: bool,
    texture: Option<egui::TextureHandle>,
    gradient: Box<dyn colorgrad::Gradient + Send + Sync>,

    // 1D curve parameters
    curve_kind: BesselKind,
    curve_scaling: Scaling,
    curve_x_range: [f64; 2],
    curve_samples: usize,
    enabled_orders: Vec<(f64, bool)>, // list of (order, is_enabled)
    custom_order: f64,
    custom_order_enabled: bool,
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

    pub fn eval_bessel(&self, x: f64, y: f64) -> f64 {
        let (order, z) = match self.domain_mode {
            DomainMode::ComplexPlane => (self.order, Complex64::new(x, y)),
            DomainMode::RealVsOrder => (y, Complex64::new(x, 0.0)),
        };

        let result = eval_bessel_single(self.kind, order, z, self.scaling);

        match result {
            Ok(val) => {
                let v = match self.component {
                    Component::Magnitude => val.norm(),
                    Component::LogMagnitude => (val.norm() + 1e-12).ln(),
                    Component::Phase => val.arg(),
                    Component::Real => val.re,
                    Component::Imag => val.im,
                };
                if v.is_finite() { v } else { f64::NAN }
            }
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

    pub fn compute_raw_grid(&self) -> Vec<f64> {
        let res = self.resolution;
        let (x_min, x_max) = (self.x_range[0], self.x_range[1]);
        let (y_min, y_max) = (self.y_range[0], self.y_range[1]);

        (0..res)
            .into_par_iter()
            .flat_map(|row| {
                let v = y_max - (row as f64 / (res - 1) as f64) * (y_max - y_min);
                (0..res)
                    .map(|col| {
                        let u = x_min + (col as f64 / (res - 1) as f64) * (x_max - x_min);
                        self.eval_bessel(u, v)
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    fn update_texture(&mut self, ctx: &egui::Context) {
        let res = self.resolution;
        let raw_values = self.compute_raw_grid();

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

        let mut pixels = Vec::with_capacity(res * res * 4);
        for &val in &raw_values {
            if val.is_nan() {
                pixels.extend_from_slice(&[30, 30, 30, 255]);
            } else {
                let t = ((val - val_min) / range_span).clamp(0.0, 1.0);
                let rgba = self.gradient.at(t as f32).to_rgba8();
                pixels.extend_from_slice(&rgba);
            }
        }

        let image = egui::ColorImage::from_rgba_unmultiplied([res, res], &pixels);
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
            let n_pts = self.curve_samples.max(50);
            let (x_min, x_max) = (self.curve_x_range[0], self.curve_x_range[1]);

            let mut active_orders = Vec::new();
            for &(order, enabled) in &self.enabled_orders {
                if enabled {
                    active_orders.push(order);
                }
            }
            if self.custom_order_enabled {
                active_orders.push(self.custom_order);
            }

            for order in active_orders {
                let points: PlotPoints = (0..n_pts)
                    .filter_map(|i| {
                        let t = i as f64 / (n_pts - 1) as f64;
                        let x = x_min + t * (x_max - x_min);
                        Self::eval_1d_bessel(self.curve_kind, order, x, self.curve_scaling)
                            .map(|y| [x, y])
                    })
                    .collect();

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

                plot_ui.line(Line::new(label, points));
            }
        });
    }
}

impl eframe::App for BesselViewer {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();

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
                        ui.label("Orders to Overlay:");
                        ui.horizontal_wrapped(|ui| {
                            for (order, enabled) in &mut self.enabled_orders {
                                ui.checkbox(enabled, format!("v = {:.0}", order));
                            }
                        });

                        ui.horizontal(|ui| {
                            ui.checkbox(&mut self.custom_order_enabled, "Custom v:");
                            if self.custom_order_enabled {
                                ui.add(egui::DragValue::new(&mut self.custom_order).speed(0.05));
                            }
                        });

                        ui.add_space(4.0);
                        ui.add(
                            egui::Slider::new(&mut self.curve_x_range[0], -50.0..=10.0)
                                .text("Curve X Min"),
                        );
                        ui.add(
                            egui::Slider::new(&mut self.curve_x_range[1], 1.0..=100.0)
                                .text("Curve X Max"),
                        );
                        ui.add(
                            egui::Slider::new(&mut self.curve_samples, 100..=2000).text("Samples"),
                        );

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
                    "2D Heatmap Surface: {}{}",
                    self.kind.label(),
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
                ui.heading(format!(
                    "1D Overlaid Curves: {}{}",
                    self.curve_kind.label(),
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
                        ui.heading(format!("2D Heatmap: {}{}", self.kind.label(), scaled_str));
                        self.show_heatmap_plot(ui);
                    });
                    columns[1].vertical(|ui| {
                        let scaled_str = if self.curve_scaling == Scaling::Scaled {
                            " [Scaled]"
                        } else {
                            ""
                        };
                        ui.heading(format!(
                            "1D Curves: {}{}",
                            self.curve_kind.label(),
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

    #[test]
    fn test_bessel_j_origin() {
        let viewer = BesselViewer {
            view_tab: ViewTab::Heatmap2D,
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
            curve_samples: 100,
            enabled_orders: vec![(0.0, true)],
            custom_order: 0.5,
            custom_order_enabled: false,
        };
        // J_0(0) = 1.0
        let val = viewer.eval_bessel(0.0, 0.0);
        assert!((val - 1.0).abs() < 1e-6, "J_0(0) should be 1.0, got {val}");
    }

    #[test]
    fn test_compute_raw_grid() {
        let viewer = BesselViewer {
            view_tab: ViewTab::Heatmap2D,
            kind: BesselKind::J,
            scaling: Scaling::Unscaled,
            domain_mode: DomainMode::ComplexPlane,
            component: Component::Magnitude,
            colormap: ColormapPreset::Turbo,
            order: 1.0,
            x_range: [-2.0, 2.0],
            y_range: [-2.0, 2.0],
            val_range: [0.0, 1.0],
            auto_range: false,
            resolution: 20,
            dirty: false,
            texture: None,
            gradient: build_gradient(ColormapPreset::Turbo),
            curve_kind: BesselKind::J,
            curve_scaling: Scaling::Unscaled,
            curve_x_range: [0.0, 20.0],
            curve_samples: 100,
            enabled_orders: vec![(0.0, true)],
            custom_order: 0.5,
            custom_order_enabled: false,
        };
        let grid = viewer.compute_raw_grid();
        assert_eq!(grid.len(), 400);
        assert!(grid.iter().any(|&v| v.is_finite()));
    }

    #[test]
    fn test_compute_raw_grid_scaled() {
        let viewer = BesselViewer {
            view_tab: ViewTab::Heatmap2D,
            kind: BesselKind::I,
            scaling: Scaling::Scaled,
            domain_mode: DomainMode::ComplexPlane,
            component: Component::Magnitude,
            colormap: ColormapPreset::Turbo,
            order: 0.0,
            x_range: [-10.0, 10.0],
            y_range: [-10.0, 10.0],
            val_range: [0.0, 1.0],
            auto_range: false,
            resolution: 20,
            dirty: false,
            texture: None,
            gradient: build_gradient(ColormapPreset::Turbo),
            curve_kind: BesselKind::I,
            curve_scaling: Scaling::Scaled,
            curve_x_range: [0.0, 20.0],
            curve_samples: 100,
            enabled_orders: vec![(0.0, true)],
            custom_order: 0.5,
            custom_order_enabled: false,
        };
        let grid = viewer.compute_raw_grid();
        assert_eq!(grid.len(), 400);
        assert!(grid.iter().all(|&v| v.is_finite() && !v.is_nan()));
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
        // Test I_0(2.0): scaled value should equal unscaled * exp(-2.0)
        let unscaled_i =
            BesselViewer::eval_1d_bessel(BesselKind::I, 0.0, 2.0, Scaling::Unscaled).unwrap();
        let scaled_i =
            BesselViewer::eval_1d_bessel(BesselKind::I, 0.0, 2.0, Scaling::Scaled).unwrap();
        let expected_scaled_i = unscaled_i * (-2.0_f64).exp();
        assert!(
            (scaled_i - expected_scaled_i).abs() < 1e-12,
            "Scaled I_0(2.0) mismatch: {scaled_i} vs expected {expected_scaled_i}"
        );

        // Test K_0(2.0): scaled value should equal unscaled * exp(2.0)
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
}
