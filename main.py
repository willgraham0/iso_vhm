import bisect
import os
import re

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from wx import (
    App,
    BoxSizer,
    Button,
    FileDialog,
    Frame,
    Panel,
    Slider,
    StaticText,
    SpinCtrl,
    RadioBox,
    ALIGN_CENTRE_VERTICAL,
    ALIGN_LEFT,
    DEFAULT_FRAME_STYLE,
    EVT_BUTTON,
    EVT_RADIOBOX,
    EVT_SPINCTRL,
    EXPAND,
    GROW,
    HORIZONTAL,
    LEFT,
    RESIZE_BORDER,
    SL_AUTOTICKS,
    SL_LABELS,
    TOP,
    VERTICAL,
    ID_OK,
    FD_OPEN,
)
from wx.lib.agw.floatspin import FloatSpin


np.seterr(invalid="ignore", divide="ignore")


class VHMEnvelope(Frame):
    gamma_rvh = 1.10

    def __init__(self, parent):
        super().__init__(
            parent, -1, "ISO V-H-M Envelope ", style=DEFAULT_FRAME_STYLE ^ RESIZE_BORDER
        )

        # Frame items.

        panel = Panel(self)
        panel.SetBackgroundColour("white")

        # Panel items.

        # Toggle fields box.
        fields = BoxSizer(HORIZONTAL)
        self.foundation = RadioBox(panel, -1, label="Foundation", choices=["Sand", "Clay"])  # Sand or Clay.
        self.unit = RadioBox(panel, -1, label="Unit", choices=["Dimensioned", "Dimensionless"])  # Units.
        self.suction = RadioBox(panel, -1, label="Suction", choices=["Available", "Unavailable"])  # Suction.
        self.alpha = RadioBox(panel, -1, label="Adhesion", choices=[u"\u0251>=0.5", u"\u0251<0.5"])  # Adhesion.
        self.factored_vh = RadioBox(panel, -1, label="Factored V-H", choices=["Off", "On"])  # Factored VH.
        fields.Add(self.foundation)
        fields.Add(self.unit)
        fields.Add(self.suction)
        fields.Add(self.alpha)
        fields.Add(self.factored_vh)

        # Qv, Vertical Force Capacity box.
        vertical_capacity = BoxSizer(HORIZONTAL)
        qv_text = StaticText(panel, -1, "Qv:")
        self.qv_value = SpinCtrl(panel, -1, min=0, max=100000, initial=10000, size=(150, -1))
        vertical_capacity.Add(qv_text, flag=ALIGN_CENTRE_VERTICAL)
        vertical_capacity.Add(self.qv_value, flag=ALIGN_CENTRE_VERTICAL)

        # Qh,Horizontal Force Capacity box.
        horizontal_capacity = BoxSizer(HORIZONTAL)
        qh_text = StaticText(panel, -1, "Qh:")
        self.qh_value = SpinCtrl(panel, -1, min=0, max=100000, initial=1000, size=(150, -1))
        horizontal_capacity.Add(qh_text, flag=ALIGN_CENTRE_VERTICAL)
        horizontal_capacity.Add(self.qh_value, flag=ALIGN_CENTRE_VERTICAL)

        # Qm, Moment Capacity box.
        moment_capacity = BoxSizer(HORIZONTAL)
        qm_text = StaticText(panel, -1, "Qm:")
        self.qm_value = SpinCtrl(panel, -1, min=0, max=100000, initial=20000, size=(150, -1))
        moment_capacity.Add(qm_text, flag=ALIGN_CENTRE_VERTICAL)
        moment_capacity.Add(self.qm_value, flag=ALIGN_CENTRE_VERTICAL)

        # Contact Diameter box.
        contact_diameter = BoxSizer(HORIZONTAL)
        diameter_text = StaticText(panel, -1, "Bmax/B:")
        self.diameter_value = FloatSpin(
            panel,
            -1,
            value=0.5,
            min_val=1.0,
            max_val=2.0,
            increment=0.1,
            digits=2,
            size=(150, -1),
            style=SL_AUTOTICKS | SL_LABELS,
        )
        contact_diameter.Add(diameter_text, flag=ALIGN_CENTRE_VERTICAL)
        contact_diameter.Add(self.diameter_value, flag=ALIGN_CENTRE_VERTICAL)

        # Depth parameter box.
        depth_parameter = BoxSizer(HORIZONTAL)
        a_text = StaticText(panel, -1, "a:")
        self.a_value = FloatSpin(
            panel,
            -1,
            value=0.0,
            min_val=0.0,
            max_val=1.0,
            increment=0.1,
            digits=2,
            size=(150, -1),
            style=SL_AUTOTICKS | SL_LABELS,
        )
        depth_parameter.Add(a_text, flag=ALIGN_CENTRE_VERTICAL)
        depth_parameter.Add(self.a_value, flag=ALIGN_CENTRE_VERTICAL)

        # Adhesion box.
        adhesion = BoxSizer(HORIZONTAL)
        alpha_text = StaticText(panel, -1, u"\u0251:")
        self.alpha_value = FloatSpin(
            panel,
            -1,
            value=0.5,
            min_val=0.5,
            max_val=1.0,
            increment=0.1,
            digits=2,
            size=(150, -1),
            style=SL_AUTOTICKS | SL_LABELS,
        )
        adhesion.Add(alpha_text, flag=ALIGN_CENTRE_VERTICAL)
        adhesion.Add(self.alpha_value, flag=ALIGN_CENTRE_VERTICAL)

        # 2-d plot slice defaults.

        # Fv, Vertical Force box.
        vertical_force = BoxSizer(HORIZONTAL)
        fv_position_text = StaticText(panel, -1, "Fv:")
        self.fv_position_value = Slider(
            panel,
            -1,
            value=self.qv_value.GetValue() // 2,
            minValue=0,
            maxValue=self.qv_value.GetValue(),
            size=(150, -1),
            style=SL_LABELS,
        )
        vertical_force.Add(fv_position_text, flag=ALIGN_CENTRE_VERTICAL | ALIGN_LEFT)
        vertical_force.Add(self.fv_position_value, flag=ALIGN_CENTRE_VERTICAL)

        # Fh, Horizontal Force box.
        horizontal_force = BoxSizer(HORIZONTAL)
        fh_position_text = StaticText(panel, -1, "Fh:")
        self.fh_value_default = 0
        self.fh_position_value = Slider(
            panel,
            -1,
            value=self.fh_value_default,
            minValue=0,
            maxValue=self.qh_value.GetValue(),
            size=(150, -1),
            style=SL_LABELS,
        )
        horizontal_force.Add(fh_position_text, flag=ALIGN_CENTRE_VERTICAL | ALIGN_LEFT)
        horizontal_force.Add(self.fh_position_value, flag=ALIGN_CENTRE_VERTICAL)

        # Fm, Moment box.
        moment = BoxSizer(HORIZONTAL)
        fm_position_text = StaticText(panel, -1, "Fm:")
        self.fm_value_default = 0
        self.fm_position_value = Slider(
            panel,
            -1,
            value=self.fm_value_default,
            minValue=0,
            maxValue=self.qm_value.GetValue(),
            size=(150, -1),
            style=SL_LABELS,
        )
        moment.Add(fm_position_text, flag=ALIGN_CENTRE_VERTICAL | ALIGN_LEFT)
        moment.Add(self.fm_position_value, flag=ALIGN_CENTRE_VERTICAL)

        # Action buttons box.
        actions = BoxSizer(HORIZONTAL)
        self.draw_button = Button(panel, -1, "Draw Plots")
        self.clear_button = Button(panel, -1, "Clear Plots and Reacts")
        self.import_button = Button(panel, -1, "Import ISO Reactions")
        self.clear_reacts_button = Button(panel, -1, "Clear Reactions")
        self.save_button = Button(panel, -1, "Save Selection")
        self.clear_save_button = Button(panel, -1, "Clear Selection")
        self.plot_selection_button = Button(panel, -1, "Plot Selection")
        actions.Add(self.draw_button, flag=ALIGN_CENTRE_VERTICAL)
        actions.Add(self.clear_button, flag=ALIGN_CENTRE_VERTICAL)
        actions.Add(self.import_button, flag=ALIGN_CENTRE_VERTICAL)
        actions.Add(self.clear_reacts_button, flag=ALIGN_CENTRE_VERTICAL)
        actions.Add(self.save_button, flag=ALIGN_CENTRE_VERTICAL)
        actions.Add(self.clear_save_button, flag=ALIGN_CENTRE_VERTICAL)
        actions.Add(self.plot_selection_button, flag=ALIGN_CENTRE_VERTICAL)

        # Wbfo box.
        wbfo = BoxSizer(HORIZONTAL)
        wbfo_text = StaticText(panel, -1, "Wbfo:")
        self.wbfo_value = SpinCtrl(
            panel, -1, min=0, max=10000, initial=200, size=(150, -1)
        )
        wbfo.Add(wbfo_text, flag=ALIGN_CENTRE_VERTICAL)
        wbfo.Add(self.wbfo_value, flag=ALIGN_CENTRE_VERTICAL)

        # Diameter box.
        factored_vh = BoxSizer(HORIZONTAL)
        self.bs_value = SpinCtrl(panel, -1, min=0, max=1000, initial=20, size=(150, -1))
        self.bs_text = StaticText(panel, -1, "Bs:")
        factored_vh.Add(self.bs_text, flag=ALIGN_CENTRE_VERTICAL)
        factored_vh.Add(self.bs_value, flag=ALIGN_CENTRE_VERTICAL)

        # Selection attributes

        self.colours = ["b", "g", "r"]
        self.moments = []
        self.horizontals = []
        self.verticals = []
        self.norm_moments = []
        self.norm_horizontals = []
        self.norm_verticals = []

        self.selection_qv = []
        self.selection_qh = []
        self.selection_qm = []
        self.selection_moments = []
        self.selection_horizontals = []
        self.selection_verticals = []
        self.selection_fh_ax1 = []
        self.selection_fv_ax1 = []
        self.selection_fh_ax2 = []
        self.selection_fm_ax2 = []
        self.selection_fv_ax3 = []
        self.selection_fm_ax3 = []
        self.selection_fh_ax4 = []
        self.selection_fv_ax4 = []
        self.selection_fm_ax4 = []
        self.selection_labels = []

        # Graph attributes

        self.dpi = 100
        self.fontSize = 8

        self.gs = gridspec.GridSpec(3, 2)
        self.fig = plt.figure(facecolor="white", dpi=self.dpi)
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)
        self.handles = None
        self.labels = None

        self.ax1 = self.fig.add_subplot(self.gs[0, 0])
        self.ax2 = self.fig.add_subplot(self.gs[1, 0])
        self.ax3 = self.fig.add_subplot(self.gs[2, 0])
        self.ax4 = self.fig.add_subplot(self.gs[0:2, 1], projection="3d")

        # Sizers

        capacities = BoxSizer(HORIZONTAL)
        capacities.Add(vertical_capacity)
        capacities.Add(horizontal_capacity)
        capacities.Add(moment_capacity)

        ratios = BoxSizer(HORIZONTAL)
        ratios.Add(contact_diameter)
        ratios.Add(depth_parameter)
        ratios.Add(adhesion)

        parameters = BoxSizer(HORIZONTAL)
        parameters.Add(capacities, flag=ALIGN_CENTRE_VERTICAL)
        parameters.Add(ratios, flag=ALIGN_CENTRE_VERTICAL)

        forces = BoxSizer(HORIZONTAL)
        forces.Add(vertical_force, flag=ALIGN_CENTRE_VERTICAL)
        forces.Add(horizontal_force, flag=ALIGN_CENTRE_VERTICAL)
        forces.Add(moment, flag=ALIGN_CENTRE_VERTICAL)
        forces.Add(wbfo, flag=ALIGN_CENTRE_VERTICAL)
        forces.Add(factored_vh, flag=ALIGN_CENTRE_VERTICAL)

        view = BoxSizer(VERTICAL)
        view.Add(fields)
        view.Add(parameters)
        view.Add(forces)
        view.Add(actions)
        view.Add(self.canvas, 1, flag=LEFT | TOP | GROW)
        view.Add(self.toolbar, flag=EXPAND)

        panel.SetSizer(view)
        view.Fit(self)

        # Events

        self.Bind(EVT_RADIOBOX, self.on_foundation_box, self.foundation)

        self.Bind(EVT_SPINCTRL, self.on_qv_spin, self.qv_value)
        self.Bind(EVT_SPINCTRL, self.on_qh_spin, self.qh_value)
        self.Bind(EVT_SPINCTRL, self.on_qm_spin, self.qm_value)

        self.Bind(EVT_BUTTON, self.draw_graphs, self.draw_button)
        self.Bind(EVT_BUTTON, self.clear_graphs, self.clear_button)
        self.Bind(EVT_BUTTON, self.import_reacts, self.import_button)
        self.Bind(EVT_BUTTON, self.save_selection, self.save_button)
        self.Bind(EVT_BUTTON, self.clear_selection, self.clear_save_button)
        self.Bind(EVT_BUTTON, self.plot_selection, self.plot_selection_button)
        self.Bind(EVT_BUTTON, self.clear_reacts_with_event, self.clear_reacts_button)

        self.Bind(EVT_RADIOBOX, self.on_factored_vh_box, self.factored_vh)

        # Initial attribute statuses

        self.a_value.SetValue(0.0)
        self.a_value.Disable()
        self.alpha_value.Disable()
        self.alpha.Disable()
        self.suction.Disable()
        self.selection_counter = 0
        self.first_plot_flag = 0

        self.wbfo_value.Disable()
        self.bs_value.Disable()

        # Initially executed methods

        self.add_graph_labels()

    # Functions

    def get_index(self, list_data, x, flag):
        i = bisect.bisect_left(list_data, x, lo=0, hi=len(list_data))
        if x == list_data[i]:
            return i
        elif x != list_data[i] and flag == "left":
            return i - 1
        elif x != list_data[i] and flag == "right":
            return i

    def lookup_table(self, table, value_a, value_alpha):
        a_values = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
        alpha_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        i = self.get_index(a_values, value_a, "left")
        as_lower = table[:, i]
        j = self.get_index(a_values, value_a, "right")
        as_upper = table[:, j]

        ii = self.get_index(alpha_values, value_alpha, "left")
        jj = self.get_index(alpha_values, value_alpha, "right")

        if i == j and ii == jj:
            value = as_lower[ii]
            return value

        if i == j and ii != jj:
            lower = as_lower[ii]
            upper = as_lower[jj]
            value = (value_alpha - alpha_values[ii]) * (upper - lower) / (
                alpha_values[jj] - alpha_values[ii]
            ) + lower
            return value

        if i != j and ii == jj:
            lower = as_lower[ii]
            upper = as_upper[ii]
            value = (value_a - a_values[i]) / (a_values[j] - a_values[i]) * (
                upper - lower
            ) + lower
            return value

        if i != j and ii != jj:
            alphas_lower = (value_a - a_values[i]) * (as_upper[ii] - as_lower[ii]) / (
                a_values[j] - a_values[i]
            ) + as_lower[ii]
            alphas_upper = (value_a - a_values[i]) * (as_upper[jj] - as_lower[jj]) / (
                a_values[j] - a_values[i]
            ) + as_lower[jj]
            value = (value_alpha - alpha_values[ii]) * (alphas_upper - alphas_lower) / (
                alpha_values[jj] - alpha_values[ii]
            ) + alphas_lower
            return value

    def make_factorisation(self, vertical_values, horizontal_values, wbfo, bs):
        verticals = vertical_values - (wbfo - bs)
        vector_lengths = np.sqrt(verticals ** 2 + horizontal_values ** 2)
        factored_vector_lengths = (1 / self.gamma_rvh) * vector_lengths

        tan_vector_lengths = verticals / vector_lengths
        cos_vector_lengths = horizontal_values / vector_lengths

        factored_verticals = factored_vector_lengths * tan_vector_lengths + (wbfo - bs)
        factored_horizontals = factored_vector_lengths * cos_vector_lengths
        return factored_verticals, factored_horizontals

    # Methods

    def clear_axes(self):
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.ax4.clear()

    def plot_selection(self, event):
        self.clear_axes()

        for i in range(len(self.selection_fh_ax1)):
            qv = self.selection_qv[i]
            qh = self.selection_qh[i]
            qm = self.selection_qm[i]

            moments = np.asarray(self.selection_moments[i])
            horizontals = np.asarray(self.selection_horizontals[i])
            verticals = np.asarray(self.selection_verticals[i])

            fh_ax1 = np.asarray(self.selection_fh_ax1[i])
            fv_ax1 = np.asarray(self.selection_fv_ax1[i])
            fh_ax2 = np.asarray(self.selection_fh_ax2[i])
            fm_ax2 = np.asarray(self.selection_fm_ax2[i])
            fv_ax3 = np.asarray(self.selection_fv_ax3[i])
            fm_ax3 = np.asarray(self.selection_fm_ax3[i])
            fh_ax4 = np.asarray(self.selection_fh_ax4[i])
            fv_ax4 = np.asarray(self.selection_fv_ax4[i])
            fm_ax4 = np.asarray(self.selection_fm_ax4[i])
            labels = self.selection_labels[i]

            if self.unit.GetSelection() == 0:
                self.ax1.plot(fh_ax1, fv_ax1, color=self.colours[i], label=labels)
                self.ax2.plot(fh_ax2, fm_ax2, color=self.colours[i])
                self.ax3.plot(fv_ax3, fm_ax3, color=self.colours[i])
                xs, ys = np.meshgrid(fv_ax4, fh_ax4)
                self.ax4.plot_surface(
                    xs, ys, fm_ax4, linewidth=0.1, alpha=0.3, color=self.colours[i]
                )

                self.ax1.scatter(
                    horizontals, verticals, color=self.colours[i], s=5, edgecolor="none"
                )
                self.ax2.scatter(
                    horizontals, moments, color=self.colours[i], s=5, edgecolor="none"
                )
                self.ax3.scatter(
                    verticals, moments, color=self.colours[i], s=5, edgecolor="none"
                )
                self.ax4.scatter(
                    verticals,
                    horizontals,
                    moments,
                    color=self.colours[i],
                    s=3,
                    linewidth="0",
                )

            if self.unit.GetSelection() == 1:
                self.ax1.plot(
                    fh_ax1 / qh, fv_ax1 / qv, color=self.colours[i], label=labels
                )
                self.ax2.plot(fh_ax2 / qh, fm_ax2 / qm, color=self.colours[i])
                self.ax3.plot(fv_ax3 / qv, fm_ax3 / qm, color=self.colours[i])
                xs, ys = np.meshgrid(fv_ax4 / qv, fh_ax4 / qh)
                self.ax4.plot_surface(
                    xs, ys, fm_ax4 / qm, linewidth=0.1, alpha=0.3, color=self.colours[i]
                )

                self.ax1.scatter(
                    horizontals / qh,
                    verticals / qv,
                    color=self.colours[i],
                    s=5,
                    edgecolor="none",
                )
                self.ax2.scatter(
                    horizontals / qh,
                    moments / qm,
                    color=self.colours[i],
                    s=5,
                    edgecolor="none",
                )
                self.ax3.scatter(
                    verticals / qv,
                    moments / qm,
                    color=self.colours[i],
                    s=5,
                    edgecolor="none",
                )
                self.ax4.scatter(
                    verticals / qv,
                    horizontals / qh,
                    moments / qm,
                    color=self.colours[i],
                    s=3,
                    linewidth="0",
                )

        self.handles, self.labels = self.ax1.get_legend_handles_labels()
        self.ax4.legend(
            self.handles,
            self.labels,
            bbox_to_anchor=(0, -0.1),
            loc="upper left",
            prop={"size": self.fontSize},
        )
        self.add_graph_labels()
        self.canvas.draw()

    def clear_selection(self, event):
        self.selection_counter = 0
        self.save_button.Enable()
        self.selection_qv = []
        self.selection_qh = []
        self.selection_qm = []
        self.selection_moments = []
        self.selection_horizontals = []
        self.selection_verticals = []
        self.selection_fh_ax1 = []
        self.selection_fv_ax1 = []
        self.selection_fh_ax2 = []
        self.selection_fm_ax2 = []
        self.selection_fv_ax3 = []
        self.selection_fm_ax3 = []
        self.selection_fh_ax4 = []
        self.selection_fv_ax4 = []
        self.selection_fm_ax4 = []
        self.selection_labels = []

    def clear_reacts_with_event(self, event):
        self.clear_reacts()

    def save_selection(self, event):
        if self.selection_counter <= 2:
            if len(self.moments) == 0:
                self.selection_moments.append([])
                self.selection_horizontals.append([])
                self.selection_verticals.append([])
            else:
                self.selection_moments.append(self.moments)
                self.selection_horizontals.append(self.horizontals)
                self.selection_verticals.append(self.verticals)

            self.selection_qv.append(self.qv_value.GetValue())
            self.selection_qh.append(self.qh_value.GetValue())
            self.selection_qm.append(self.qm_value.GetValue())

            self.selection_fh_ax1.append(self.fh_ax1)
            self.selection_fv_ax1.append(self.fv_ax1)
            self.selection_fh_ax2.append(self.fh_ax2)
            self.selection_fm_ax2.append(self.fm_ax2)
            self.selection_fv_ax3.append(self.fv_ax3)
            self.selection_fm_ax3.append(self.fm_ax3)
            self.selection_fh_ax4.append(self.fh_ax4)
            self.selection_fv_ax4.append(self.fv_ax4)
            self.selection_fm_ax4.append(self.fm_ax4)
            self.selection_labels.append(self.label)

            self.selection_counter = self.selection_counter + 1
            if self.selection_counter > 2:
                self.save_button.Disable()

    def on_foundation_box(self, event):
        if self.foundation.GetSelection() == 0:
            self.suction.Disable()
            self.alpha.Disable()
            self.a_value.SetValue(0.00)
            self.alpha_value.SetValue(0.5)
            self.a_value.Disable()
            self.alpha_value.Disable()
            self.diameter_value.SetValue(1.00)
            self.diameter_value.Enable()
        else:
            self.diameter_value.SetValue(1.00)
            self.diameter_value.Disable()
            self.a_value.Enable()
            self.alpha_value.Enable()
            self.alpha.Enable()
            self.suction.Enable()

    def on_factored_vh_box(self, event):
        if self.factored_vh.GetSelection() == 0:
            self.wbfo_value.Disable()
            self.bs_value.Disable()
        else:
            self.wbfo_value.Enable()
            self.bs_value.Enable()

    def on_qv_spin(self, event):
        self.fv_position_value.SetMax(self.qv_value.GetValue())
        self.fv_position_value.SetValue(self.qv_value.GetValue() / 2)

    def on_qh_spin(self, event):
        self.fh_position_value.SetMax(self.qh_value.GetValue())
        self.fh_position_value.SetValue(self.fh_value_default)

    def on_qm_spin(self, event):
        self.fm_position_value.SetMax(self.qm_value.GetValue())
        self.fm_position_value.SetValue(self.fm_value_default)

    def import_reacts(self, event):

        dialog = FileDialog(
            None, "Choose a .JAR file", defaultDir=os.getcwd(), style=FD_OPEN
        )
        if dialog.ShowModal() == ID_OK:
            path = dialog.GetPath()
            inFile = open(path, "r")
            data = inFile.read()
            inFile.close()
            reactions = re.findall(
                r"\d{1,3}\.\d\s*[a-d]\s*\d\s*[\-0]\.\d{4}E[\+\-]\d{2}\s*\-\d*\.\s*\d*\.\s*([\-0]\.\d{4}E[\+\-]\d{2})\s*\-(\d*\.)\s*(\d*\.)\s*\d*\.\s*\d\.\d{2}\s*\d\s*[\-0]\.\d{4}E[\+\-]\d{2}\s*\-\d*\.\s*\d*\.\s*([\-0]\.\d{4}E[\+\-]\d{2})\s*\-(\d*\.)\s*(\d*\.)\s*\d*\.\s*\d\.\d{2}\s*\d\s*[\-0]\.\d{4}E[\+\-]\d{2}\s*\-\d*\.\s*\d*\.\s*([\-0]\.\d{4}E[\+\-]\d{2})\s*\-(\d*\.)\s*(\d*\.)\s*\d*\.\s*\d\.\d{2}",
                data,
            )

            self.moments = np.asarray(
                [[float(i[0]), float(i[3]), float(i[6])] for i in reactions]
            ).flatten()
            self.horizontals = np.asarray(
                [[float(i[1]), float(i[4]), float(i[7])] for i in reactions]
            ).flatten()
            self.verticals = np.asarray(
                [[float(i[2]), float(i[5]), float(i[8])] for i in reactions]
            ).flatten()

    def normalise_reacts(self):

        self.norm_moments = np.array(self.moments) / float(self.qm_value.GetValue())
        self.norm_horizontals = np.array(self.horizontals) / float(
            self.qh_value.GetValue()
        )
        self.norm_verticals = np.array(self.verticals) / float(self.qv_value.GetValue())

    def draw_reacts(self):

        if all([self.horizontals, self.verticals, self.moments]):
            self.normalise_reacts()

            if self.unit.GetSelection() == 0:
                self.ax1.scatter(
                    self.horizontals, self.verticals, s=5, edgecolor="none"
                )
                self.ax2.scatter(self.horizontals, self.moments, s=5, edgecolor="none")
                self.ax3.scatter(self.verticals, self.moments, s=5, edgecolor="none")
                self.ax4.scatter(
                    self.verticals, self.horizontals, self.moments, s=5, linewidth="0"
                )
                self.canvas.draw()
            else:
                self.ax1.scatter(
                    self.norm_horizontals, self.norm_verticals, s=5, edgecolor="none"
                )
                self.ax2.scatter(
                    self.norm_horizontals, self.norm_moments, s=5, edgecolor="none"
                )
                self.ax3.scatter(
                    self.norm_verticals, self.norm_moments, s=5, edgecolor="none"
                )
                self.ax4.scatter(
                    self.norm_verticals,
                    self.norm_horizontals,
                    self.norm_moments,
                    s=5,
                    linewidth="0",
                )
                self.canvas.draw()

            self.add_graph_labels()

    def add_graph_labels(self):
        self.ax1.grid(True)
        self.ax2.grid(True)
        self.ax3.grid(True)
        self.ax4.grid(False)

        self.ax1.tick_params(axis="both", which="major", labelsize=self.fontSize)
        self.ax2.tick_params(axis="both", which="major", labelsize=self.fontSize)
        self.ax3.tick_params(axis="both", which="major", labelsize=self.fontSize)
        self.ax4.tick_params(axis="both", which="major", labelsize=self.fontSize)

        if self.unit.GetSelection() == 0:
            self.ax1.set_xlabel(r"$\mathregular{F_H}$ (tonnes)", fontsize=self.fontSize)
            self.ax2.set_xlabel(r"$\mathregular{F_H}$ (tonnes)", fontsize=self.fontSize)
            self.ax3.set_xlabel(r"$\mathregular{F_V}$ (tonnes)", fontsize=self.fontSize)
            self.ax4.set_xlabel(r"$\mathregular{F_V}$ (tonnes)", fontsize=self.fontSize)

            self.ax1.set_ylabel(r"$\mathregular{F_V}$ (tonnes)", fontsize=self.fontSize)
            self.ax2.set_ylabel(r"$\mathregular{F_M}$ (tonnes)", fontsize=self.fontSize)
            self.ax3.set_ylabel(r"$\mathregular{F_M}$ (tonnes)", fontsize=self.fontSize)
            self.ax4.set_ylabel(r"$\mathregular{F_H}$ (tonnes)", fontsize=self.fontSize)

            self.ax4.set_zlabel(
                r"Moment Reaction, $\mathregular{F_M}$ (tonnes)", fontsize=self.fontSize
            )
        else:
            self.ax1.set_xlabel(r"$\mathregular{F_H/Q_H}$", fontsize=self.fontSize)
            self.ax2.set_xlabel(r"$\mathregular{F_H/Q_H}$", fontsize=self.fontSize)
            self.ax3.set_xlabel(r"$\mathregular{F_V/Q_V}$", fontsize=self.fontSize)
            self.ax4.set_xlabel(r"$\mathregular{F_V/Q_V}$", fontsize=self.fontSize)

            self.ax1.set_ylabel(r"$\mathregular{F_V/Q_V}$", fontsize=self.fontSize)
            self.ax2.set_ylabel(r"$\mathregular{F_M/Q_M}$", fontsize=self.fontSize)
            self.ax3.set_ylabel(r"$\mathregular{F_M/Q_M}$", fontsize=self.fontSize)
            self.ax4.set_ylabel(r"$\mathregular{F_H/Q_H}$", fontsize=self.fontSize)

            self.ax4.set_zlabel(r"$\mathregular{F_M/Q_M}$", fontsize=self.fontSize)

        self.ax1.set_xlim(left=0)
        self.ax1.set_ylim(bottom=0)
        self.ax2.set_xlim(left=0)
        self.ax2.set_ylim(bottom=0)
        self.ax3.set_xlim(left=0)
        self.ax3.set_ylim(bottom=0)

    def clear_reacts(self):
        self.moments = []
        self.horizontals = []
        self.verticals = []
        self.norm_moments = []
        self.norm_horizontals = []
        self.norm_verticals = []

    def clear_graphs(self, event):
        self.clear_axes()
        self.clear_reacts()
        self.add_graph_labels()
        self.canvas.draw()

    def draw_graphs(self, event):
        self.clear_axes()

        fv_qvts = np.array(
            [
                [0.354, 0.334, 0.308, 0.293, 0.276, 0.238, 0.200],
                [0.387, 0.373, 0.354, 0.343, 0.331, 0.300, 0.265],
                [0.418, 0.408, 0.396, 0.388, 0.379, 0.357, 0.329],
                [0.447, 0.441, 0.433, 0.428, 0.423, 0.409, 0.390],
                [0.474, 0.471, 0.468, 0.465, 0.463, 0.457, 0.448],
                [0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500],
            ]
        )

        m_alphas = np.array(
            [
                [1.172, 1.200, 1.239, 1.264, 1.295, 1.378, 1.500],
                [0.902, 0.917, 0.937, 0.950, 0.965, 1.006, 1.067],
                [0.653, 0.661, 0.670, 0.676, 0.683, 0.701, 0.729],
                [0.422, 0.425, 0.429, 0.431, 0.434, 0.440, 0.450],
                [0.205, 0.206, 0.207, 0.207, 0.208, 0.209, 0.211],
                [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
            ]
        )

        granularity = 300
        qv = float(self.qv_value.GetValue())
        qh = float(self.qh_value.GetValue())
        qm = float(self.qm_value.GetValue())
        bm_b = float(self.diameter_value.GetValue())
        a = float(self.a_value.GetValue())
        alpha = float(self.alpha_value.GetValue())
        alpha_flag = self.alpha.GetSelection()
        fv_position = float(self.fv_position_value.GetValue())
        fh_position = float(self.fh_position_value.GetValue())
        fm_position = float(self.fm_position_value.GetValue())
        wbfo = float(self.wbfo_value.GetValue())
        bs = float(self.bs_value.GetValue())

        if self.foundation.GetSelection() == 0:
            soil_flag = "Sand"
        else:
            soil_flag = "Clay"
            suction_flag = self.suction.GetSelection()

        if self.unit.GetSelection() == 0:
            graph_flag = "Dimensioned"
        else:
            graph_flag = "Dimensionless"

        startv = 0.0
        starth = 0.0
        fv_qvt = self.lookup_table(fv_qvts, a, alpha)
        m_alpha = self.lookup_table(m_alphas, a, alpha)

        if soil_flag == "Sand":
            self.fv_ax1 = np.linspace(startv, qv, granularity)
            qmps_ax1 = qm * (bm_b) ** 3 * np.ones(granularity)
            qmpv_ax1 = 0.15 * ((qm * 0.12) / (0.075 * qh)) * self.fv_ax1
            fh_ax1_orig = np.nan_to_num(
                qh
                * np.sqrt(
                    4 * a * (self.fv_ax1 / qv) * (1 - (self.fv_ax1 / qv))
                    + 16
                    * (1 - a)
                    * (self.fv_ax1 / qv) ** 2
                    * (1 - (self.fv_ax1 / qv)) ** 2
                    - (fm_position / qm) ** 2
                )
            )
            fh_ax1_mod = np.nan_to_num(
                qh
                * np.sqrt(
                    4 * a * (self.fv_ax1 / qv) * (1 - (self.fv_ax1 / qv))
                    + 16
                    * (1 - a)
                    * (self.fv_ax1 / qv) ** 2
                    * (1 - (self.fv_ax1 / qv)) ** 2
                    - (fm_position / np.minimum(qmps_ax1, qmpv_ax1)) ** 2
                )
            )
            self.fh_ax1 = np.maximum(fh_ax1_orig, fh_ax1_mod)

            if fm_position == 0:
                self.fv_ax1_factored, self.fh_ax1_factored = self.make_factorisation(
                    self.fv_ax1, self.fh_ax1, wbfo, bs
                )

            self.fh_ax2 = np.linspace(starth, qh, granularity)
            qmps_ax2 = qm * (bm_b) ** 3 * np.ones(granularity)
            qmpv_ax2 = 0.15 * ((qm * 0.12) / (0.075 * qh)) * fv_position
            fm_ax2_orig = np.nan_to_num(
                qm
                * np.sqrt(
                    4 * a * (fv_position / qv) * (1 - (fv_position / qv))
                    + 16
                    * (1 - a)
                    * (fv_position / qv) ** 2
                    * (1 - (fv_position / qv)) ** 2
                    - (self.fh_ax2 / qh) ** 2
                )
            )
            fm_ax2_mod = np.nan_to_num(
                np.minimum(qmps_ax2, qmpv_ax2)
                * np.sqrt(
                    4 * a * (fv_position / qv) * (1 - (fv_position / qv))
                    + 16
                    * (1 - a)
                    * (fv_position / qv) ** 2
                    * (1 - (fv_position / qv)) ** 2
                    - (self.fh_ax2 / qh) ** 2
                )
            )
            self.fm_ax2 = np.maximum(fm_ax2_orig, fm_ax2_mod)

            self.fv_ax3 = np.linspace(startv, qv, granularity)
            qmps_ax3 = qm * (bm_b) ** 3 * np.ones(granularity)
            qmpv_ax3 = 0.15 * ((qm * 0.12) / (0.075 * qh)) * self.fv_ax3
            fm_ax3_orig = qm * np.sqrt(
                4 * a * (self.fv_ax3 / qv) * (1 - (self.fv_ax3 / qv))
                + 16 * (1 - a) * (self.fv_ax3 / qv) ** 2 * (1 - (self.fv_ax3 / qv)) ** 2
                - (fh_position / qh) ** 2
            )
            fm_ax3_mod = np.minimum(qmps_ax3, qmpv_ax3) * np.sqrt(
                4 * a * (self.fv_ax3 / qv) * (1 - (self.fv_ax3 / qv))
                + 16 * (1 - a) * (self.fv_ax3 / qv) ** 2 * (1 - (self.fv_ax3 / qv)) ** 2
                - (fh_position / qh) ** 2
            )
            self.fm_ax3 = np.maximum(fm_ax3_orig, fm_ax3_mod)

            self.fv_ax4 = np.linspace(startv, qv, granularity)
            self.fh_ax4 = np.linspace(starth, qh, granularity)
            qmps_ax4 = qm * (bm_b) ** 3 * np.ones(granularity)
            qmpv_ax4 = 0.15 * ((qm * 0.12) / (0.075 * qh)) * self.fv_ax4
            fm_ax4_before = []

            for i in range(granularity):
                fm_inter_orig = qm * np.sqrt(
                    4 * a * (self.fv_ax4 / qv) * (1 - self.fv_ax4 / qv)
                    + 16
                    * (1 - a)
                    * (self.fv_ax4 / qv) ** 2
                    * (1 - self.fv_ax4 / qv) ** 2
                    - (self.fh_ax4[i] / qh) ** 2
                )
                fm_inter_mod = np.minimum(qmps_ax4, qmpv_ax4) * np.sqrt(
                    4 * a * (self.fv_ax4 / qv) * (1 - self.fv_ax4 / qv)
                    + 16
                    * (1 - a)
                    * (self.fv_ax4 / qv) ** 2
                    * (1 - self.fv_ax4 / qv) ** 2
                    - (self.fh_ax4[i] / qh) ** 2
                )
                fm_inter = np.maximum(fm_inter_orig, fm_inter_mod)
                fm_ax4_before.append(fm_inter)

            self.fm_ax4 = np.nan_to_num(fm_ax4_before)
            xs, ys = np.meshgrid(self.fv_ax4, self.fh_ax4)

            self.label = (
                "$\mathregular{Q_V}$=%dt; $\mathregular{Q_H}$=%dt; $\mathregular{Q_M}$=%dt; $\mathregular{B_{max}}$/B= %.2f"
                % (qv, qh, qm, bm_b)
            )

            if graph_flag == "Dimensioned":
                self.ax1.plot(self.fh_ax1, self.fv_ax1, label=self.label)
                self.ax2.plot(self.fh_ax2, self.fm_ax2)
                self.ax3.plot(self.fv_ax3, self.fm_ax3)
                self.ax4.plot_surface(xs, ys, self.fm_ax4, linewidth=0.1, alpha=0.3)

                if self.factored_vh.GetSelection() == 1:
                    if fm_position == 0:
                        self.ax1.plot(
                            self.fh_ax1_factored, self.fv_ax1_factored, "--", color="b"
                        )
                        self.ax1.plot(0, wbfo - bs, ".", ms=5, color="b")

            elif graph_flag == "Dimensionless":
                self.ax1.plot(self.fh_ax1 / qh, self.fv_ax1 / qv, label=self.label)
                self.ax2.plot(self.fh_ax2 / qh, self.fm_ax2 / qm)
                self.ax3.plot(self.fv_ax3 / qv, self.fm_ax3 / qm)
                self.ax4.plot_surface(
                    xs / qv, ys / qh, self.fm_ax4 / qm, linewidth=0.1, alpha=0.3
                )

                if self.factored_vh.GetSelection() == 1:
                    if fm_position == 0:
                        self.ax1.plot(
                            self.fh_ax1_factored / qh,
                            self.fv_ax1_factored / qv,
                            "--",
                            color="b",
                        )
                        self.ax1.plot(0, (wbfo - bs) / qv, ".", ms=5, color="b")

        if soil_flag == "Clay":
            if alpha_flag == 1:
                self.fv_ax1 = np.linspace(startv, qv, granularity)
                self.fh_ax1 = qh * np.sqrt(
                    4 * a * (self.fv_ax1 / qv) * (1 - (self.fv_ax1 / qv))
                    + 16
                    * (1 - a)
                    * (self.fv_ax1 / qv) ** 2
                    * (1 - (self.fv_ax1 / qv)) ** 2
                    - (fm_position / qm) ** 2
                )

                if fm_position == 0:
                    self.fv_ax1_factored, self.fh_ax1_factored = self.make_factorisation(
                        self.fv_ax1, self.fh_ax1, wbfo, bs
                    )

                self.fh_ax2 = np.linspace(starth, qh, granularity)
                self.fm_ax2 = qm * np.sqrt(
                    4 * a * (fv_position / qv) * (1 - (fv_position / qv))
                    + 16
                    * (1 - a)
                    * (fv_position / qv) ** 2
                    * (1 - (fv_position / qv)) ** 2
                    - (self.fh_ax2 / qh) ** 2
                )

                self.fv_ax3 = np.linspace(startv, qv, granularity)
                self.fm_ax3 = qm * np.sqrt(
                    4 * a * (self.fv_ax3 / qv) * (1 - (self.fv_ax3 / qv))
                    + 16
                    * (1 - a)
                    * (self.fv_ax3 / qv) ** 2
                    * (1 - (self.fv_ax3 / qv)) ** 2
                    - (fh_position / qh) ** 2
                )

                self.fv_ax4 = np.linspace(startv, qv, granularity)
                self.fh_ax4 = np.linspace(starth, qh, granularity)
                fm_ax4_before = []

                for i in range(granularity):
                    fm_inter = qm * np.sqrt(
                        4 * a * (self.fv_ax4 / qv) * (1 - self.fv_ax4 / qv)
                        + 16
                        * (1 - a)
                        * (self.fv_ax4 / qv) ** 2
                        * (1 - self.fv_ax4 / qv) ** 2
                        - (self.fh_ax4[i] / qh) ** 2
                    )
                    fm_ax4_before.append(fm_inter)
                self.fm_ax4 = np.nan_to_num(fm_ax4_before)
                xs, ys = np.meshgrid(self.fv_ax4, self.fh_ax4)

                self.label = (
                    "$\mathregular{Q_V}$=%dt; $\mathregular{Q_H}$=%dt; $\mathregular{Q_M}$=%dt; $\mathregular{\\alpha}$<0.50; a= %.2f; Suction N/A"
                    % (qv, qh, qm, a)
                )

                if graph_flag == "Dimensioned":
                    self.ax1.plot(self.fh_ax1, self.fv_ax1, label=self.label)
                    self.ax2.plot(self.fh_ax2, self.fm_ax2)
                    self.ax3.plot(self.fv_ax3, self.fm_ax3)
                    self.ax4.plot_surface(xs, ys, self.fm_ax4, linewidth=0.1, alpha=0.3)

                    if self.factored_vh.GetSelection() == 1:
                        if fm_position == 0:
                            self.ax1.plot(
                                self.fh_ax1_factored,
                                self.fv_ax1_factored,
                                "--",
                                color="b",
                            )
                            self.ax1.plot(0, wbfo - bs, ".", ms=5, color="b")

                elif graph_flag == "Dimensionless":
                    self.ax1.plot(self.fh_ax1 / qh, self.fv_ax1 / qv, label=self.label)
                    self.ax2.plot(self.fh_ax2 / qh, self.fm_ax2 / qm)
                    self.ax3.plot(self.fv_ax3 / qv, self.fm_ax3 / qm)
                    self.ax4.plot_surface(
                        xs / qv, ys / qh, self.fm_ax4 / qm, linewidth=0.1, alpha=0.3
                    )

                if self.factored_vh.GetSelection() == 1:
                    if fm_position == 0:
                        self.ax1.plot(
                            self.fh_ax1_factored / qh,
                            self.fv_ax1_factored / qv,
                            "--",
                            color="b",
                        )
                        self.ax1.plot(0, (wbfo - bs) / qv, ".", ms=5, color="b")

            if alpha_flag == 0:
                if suction_flag == 0:
                    fv_ax1_part1 = qv * np.linspace(0.0, fv_qvt, granularity)
                    fv_ax1_part2 = np.linspace(qv * fv_qvt, qv, granularity)
                    self.fv_ax1 = np.append(fv_ax1_part1, fv_ax1_part2)
                    fh_ax1_part1 = qh * np.sqrt(
                        (alpha + m_alpha * (fv_ax1_part1 / qv)) ** 2
                        - (fm_position / qm) ** 2
                    )
                    fh_ax1_part2 = qh * np.sqrt(
                        4 * a * (fv_ax1_part2 / qv) * (1 - (fv_ax1_part2 / qv))
                        + 16
                        * (1 - a)
                        * (fv_ax1_part2 / qv) ** 2
                        * (1 - (fv_ax1_part2 / qv)) ** 2
                        - (fm_position / qm) ** 2
                    )
                    self.fh_ax1 = np.append(fh_ax1_part1, fh_ax1_part2)

                    if fm_position == 0:
                        self.fv_ax1_factored, self.fh_ax1_factored = self.make_factorisation(
                            self.fv_ax1, self.fh_ax1, wbfo, bs
                        )

                    fh_part1_lim = qh * np.sqrt(
                        (alpha + m_alpha * fv_qvt) ** 2 - (0 / qm) ** 2
                    )
                    fh_ax2_part1 = np.linspace(0.0, fh_part1_lim, granularity)
                    fh_ax2_part2 = np.linspace(0.0, qh, granularity)
                    fm_ax2_part1 = qm * np.sqrt(
                        (alpha + m_alpha * (fv_position / qv)) ** 2
                        - (fh_ax2_part1 / qh) ** 2
                    )
                    fm_ax2_part2 = qm * np.sqrt(
                        4 * a * (fv_position / qv) * (1 - (fv_position / qv))
                        + 16
                        * (1 - a)
                        * (fv_position / qv) ** 2
                        * (1 - (fv_position / qv)) ** 2
                        - (fh_ax2_part2 / qh) ** 2
                    )

                    if fv_position <= qv * fv_qvt:
                        self.fh_ax2 = fh_ax2_part1
                        self.fm_ax2 = fm_ax2_part1
                    elif fv_position > qv * fv_qvt:
                        self.fh_ax2 = fh_ax2_part2
                        self.fm_ax2 = fm_ax2_part2

                    fv_ax3_part1 = qv * np.linspace(0.0, fv_qvt, granularity)
                    fv_ax3_part2 = np.linspace(qv * fv_qvt, qv, granularity)
                    self.fv_ax3 = np.append(fv_ax3_part1, fv_ax3_part2)
                    fm_ax3_part1 = qm * np.sqrt(
                        (alpha + m_alpha * (fv_ax3_part1 / qv)) ** 2
                        - (fh_position / qh) ** 2
                    )
                    fm_ax3_part2 = qm * np.sqrt(
                        4 * a * (fv_ax3_part2 / qv) * (1 - (fv_ax3_part2 / qv))
                        + 16
                        * (1 - a)
                        * (fv_ax3_part2 / qv) ** 2
                        * (1 - (fv_ax3_part2 / qv)) ** 2
                        - (fh_position / qh) ** 2
                    )
                    self.fm_ax3 = np.append(fm_ax3_part1, fm_ax3_part2)

                    fv_ax4_part1 = np.linspace(0.0, qv * fv_qvt, granularity // 2)
                    self.fh_ax4 = np.linspace(0.0, qh, granularity)
                    fm_ax4_part1 = []

                    for i in range(granularity):
                        fm_ax4_part1_inter = qm * np.sqrt(
                            (alpha + m_alpha * (fv_ax4_part1 / qv)) ** 2
                            - (self.fh_ax4[i] / qh) ** 2
                        )
                        fm_ax4_part1.append(fm_ax4_part1_inter)
                    fm_ax4_part1_nan = np.nan_to_num(fm_ax4_part1)

                    fv_ax4_part2 = np.linspace(qv * fv_qvt, qv, granularity // 2)
                    fm_ax4_part2 = []

                    for i in range(granularity):
                        fm_ax4_part2_inter = qm * np.sqrt(
                            4 * a * (fv_ax4_part2 / qv) * (1 - fv_ax4_part2 / qv)
                            + 16
                            * (1 - a)
                            * (fv_ax4_part2 / qv) ** 2
                            * (1 - fv_ax4_part2 / qv) ** 2
                            - (self.fh_ax4[i] / qh) ** 2
                        )
                        fm_ax4_part2.append(fm_ax4_part2_inter)
                    fm_ax4_part2_nan = np.nan_to_num(fm_ax4_part2)

                    self.fv_ax4 = np.append(fv_ax4_part1, fv_ax4_part2)
                    self.fm_ax4 = np.append(fm_ax4_part1_nan, fm_ax4_part2_nan, axis=1)
                    Xs, Ys = np.meshgrid(self.fv_ax4, self.fh_ax4)

                    self.label = (
                        "$\mathregular{Q_V}$=%dt; $\mathregular{Q_H}$=%dt; $\mathregular{Q_M}$=%dt; $\mathregular{\\alpha}$= %.2f; a= %.2f; Suction"
                        % (qv, qh, qm, alpha, a)
                    )

                    if graph_flag == "Dimensioned":
                        self.ax1.plot(self.fh_ax1, self.fv_ax1, label=self.label)
                        self.ax2.plot(self.fh_ax2, self.fm_ax2)
                        self.ax3.plot(self.fv_ax3, self.fm_ax3)
                        self.ax4.plot_surface(
                            Xs, Ys, self.fm_ax4, linewidth=0.1, alpha=0.3
                        )

                        if self.factored_vh.GetSelection() == 1:
                            if fm_position == 0:
                                self.ax1.plot(
                                    self.fh_ax1_factored,
                                    self.fv_ax1_factored,
                                    "--",
                                    color="b",
                                )
                                self.ax1.plot(0, wbfo - bs, ".", ms=5, color="b")

                    elif graph_flag == "Dimensionless":
                        self.ax1.plot(
                            self.fh_ax1 / qh, self.fv_ax1 / qv, label=self.label
                        )
                        self.ax2.plot(self.fh_ax2 / qh, self.fm_ax2 / qm)
                        self.ax3.plot(self.fv_ax3 / qv, self.fm_ax3 / qm)
                        self.ax4.plot_surface(
                            Xs / qv, Ys / qh, self.fm_ax4 / qm, linewidth=0.1, alpha=0.3
                        )

                        if self.factored_vh.GetSelection() == 1:
                            if fm_position == 0:
                                self.ax1.plot(
                                    self.fh_ax1_factored / qh,
                                    self.fv_ax1_factored / qv,
                                    "--",
                                    color="b",
                                )
                                self.ax1.plot(0, (wbfo - bs) / qv, ".", ms=5, color="b")

                if suction_flag == 1:
                    fv_ax1_part1 = np.linspace(1, qv * fv_qvt, granularity)
                    fv_ax1_part2 = np.linspace(qv * fv_qvt, qv, granularity)
                    self.fv_ax1 = np.append(fv_ax1_part1, fv_ax1_part2)
                    fh_ax1_part1 = (
                        qh
                        * (alpha + m_alpha * (fv_ax1_part1 / qv))
                        * np.sqrt(
                            1
                            - (
                                fm_position
                                / (
                                    qm
                                    * np.sqrt(
                                        16
                                        * (1 - a)
                                        * (fv_ax1_part1 / qv) ** 2
                                        * (1 - (fv_ax1_part1 / qv)) ** 2
                                        + 4
                                        * a
                                        * (fv_ax1_part1 / qv)
                                        * (1 - (fv_ax1_part1 / qv))
                                    )
                                )
                            )
                            ** 2
                        )
                    )
                    fh_ax1_part2 = qh * np.sqrt(
                        4 * a * (fv_ax1_part2 / qv) * (1 - (fv_ax1_part2 / qv))
                        + 16
                        * (1 - a)
                        * (fv_ax1_part2 / qv) ** 2
                        * (1 - (fv_ax1_part2 / qv)) ** 2
                        - (fm_position / qm) ** 2
                    )
                    self.fh_ax1 = np.append(fh_ax1_part1, fh_ax1_part2)

                    if fm_position == 0:
                        self.fv_ax1_factored, self.fh_ax1_factored = self.make_factorisation(
                            self.fv_ax1, self.fh_ax1, wbfo, bs
                        )

                    fh_part1_lim = qh * np.sqrt(
                        (alpha + m_alpha * fv_qvt) ** 2 - (0 / qm) ** 2
                    )
                    fh_ax2_part1 = np.linspace(0.0, fh_part1_lim, granularity)
                    fh_ax2_part2 = np.linspace(0.0, qh, granularity)
                    fm_ax2_part1 = (
                        qm
                        * np.sqrt(
                            16
                            * (1 - a)
                            * (fv_position / qv) ** 2
                            * (1 - (fv_position / qv)) ** 2
                            + 4 * a * (fv_position / qv) * (1 - (fv_position / qv))
                        )
                    ) * np.sqrt(
                        1
                        - (fh_ax2_part1 / ((alpha + m_alpha * (fv_position / qv)) * qh))
                        ** 2
                    )
                    fm_ax2_part2 = qm * np.sqrt(
                        4 * a * (fv_position / qv) * (1 - (fv_position / qv))
                        + 16
                        * (1 - a)
                        * (fv_position / qv) ** 2
                        * (1 - (fv_position / qv)) ** 2
                        - (fh_ax2_part2 / qh) ** 2
                    )

                    if fv_position <= qv * fv_qvt:
                        self.fh_ax2 = fh_ax2_part1
                        self.fm_ax2 = fm_ax2_part1
                    elif fv_position > qv * fv_qvt:
                        self.fh_ax2 = fh_ax2_part2
                        self.fm_ax2 = fm_ax2_part2

                    fv_ax3_part1 = qv * np.linspace(0.0, fv_qvt, granularity)
                    fv_ax3_part2 = np.linspace(qv * fv_qvt, qv, granularity)
                    self.fv_ax3 = np.append(fv_ax3_part1, fv_ax3_part2)
                    fm_ax3_part1 = (
                        qm
                        * np.sqrt(
                            16
                            * (1 - a)
                            * (fv_ax3_part1 / qv) ** 2
                            * (1 - (fv_ax3_part1 / qv)) ** 2
                            + 4 * a * (fv_ax3_part1 / qv) * (1 - (fv_ax3_part1 / qv))
                        )
                    ) * np.sqrt(
                        1
                        - (fh_position / ((alpha + m_alpha * (fv_ax3_part1 / qv)) * qh))
                        ** 2
                    )
                    fm_ax3_part2 = qm * np.sqrt(
                        4 * a * (fv_ax3_part2 / qv) * (1 - (fv_ax3_part2 / qv))
                        + 16
                        * (1 - a)
                        * (fv_ax3_part2 / qv) ** 2
                        * (1 - (fv_ax3_part2 / qv)) ** 2
                        - (fh_position / qh) ** 2
                    )
                    self.fm_ax3 = np.append(fm_ax3_part1, fm_ax3_part2)

                    fv_ax4_part1 = np.linspace(0.0, qv * fv_qvt, granularity // 2)
                    self.fh_ax4 = np.linspace(0.0, qh, granularity)
                    fm_ax4_part1 = []

                    for i in range(granularity):
                        fm_ax4_part1_inter = (
                            qm
                            * np.sqrt(
                                16
                                * (1 - a)
                                * (fv_ax4_part1 / qv) ** 2
                                * (1 - (fv_ax4_part1 / qv)) ** 2
                                + 4
                                * a
                                * (fv_ax4_part1 / qv)
                                * (1 - (fv_ax4_part1 / qv))
                            )
                        ) * np.sqrt(
                            1
                            - (
                                self.fh_ax4[i]
                                / ((alpha + m_alpha * (fv_ax4_part1 / qv)) * qh)
                            )
                            ** 2
                        )
                        fm_ax4_part1.append(fm_ax4_part1_inter)
                    fm_ax4_part1_nan = np.nan_to_num(fm_ax4_part1)

                    fv_ax4_part2 = np.linspace(qv * fv_qvt, qv, granularity // 2)
                    fm_ax4_part2 = []

                    for i in range(granularity):
                        fm_ax4_part2_inter = qm * np.sqrt(
                            4 * a * (fv_ax4_part2 / qv) * (1 - fv_ax4_part2 / qv)
                            + 16
                            * (1 - a)
                            * (fv_ax4_part2 / qv) ** 2
                            * (1 - fv_ax4_part2 / qv) ** 2
                            - (self.fh_ax4[i] / qh) ** 2
                        )
                        fm_ax4_part2.append(fm_ax4_part2_inter)
                    fm_ax4_part2_nan = np.nan_to_num(fm_ax4_part2)

                    self.fv_ax4 = np.append(fv_ax4_part1, fv_ax4_part2)
                    self.fm_ax4 = np.append(fm_ax4_part1_nan, fm_ax4_part2_nan, axis=1)
                    Xs, Ys = np.meshgrid(self.fv_ax4, self.fh_ax4)

                    self.label = (
                        "$\mathregular{Q_V}$=%dt; $\mathregular{Q_H}$=%dt; $\mathregular{Q_M}$=%dt; $\mathregular{\\alpha}$= %.2f; a= %.2f; No Suction"
                        % (qv, qh, qm, alpha, a)
                    )
                    if graph_flag == "Dimensioned":
                        self.ax1.plot(self.fh_ax1, self.fv_ax1, label=self.label)
                        self.ax2.plot(self.fh_ax2, self.fm_ax2)
                        self.ax3.plot(self.fv_ax3, self.fm_ax3)
                        self.ax4.plot_surface(
                            Xs, Ys, self.fm_ax4, linewidth=0.1, alpha=0.3
                        )

                        if self.factored_vh.GetSelection() == 1:
                            if fm_position == 0:
                                self.ax1.plot(
                                    self.fh_ax1_factored,
                                    self.fv_ax1_factored,
                                    "--",
                                    color="b",
                                )
                                self.ax1.plot(0, wbfo - bs, ".", ms=5, color="b")

                    elif graph_flag == "Dimensionless":
                        self.ax1.plot(
                            self.fh_ax1 / qh, self.fv_ax1 / qv, label=self.label
                        )
                        self.ax2.plot(self.fh_ax2 / qh, self.fm_ax2 / qm)
                        self.ax3.plot(self.fv_ax3 / qv, self.fm_ax3 / qm)
                        self.ax4.plot_surface(
                            Xs / qv, Ys / qh, self.fm_ax4 / qm, linewidth=0.1, alpha=0.3
                        )

                        if self.factored_vh.GetSelection() == 1:
                            if fm_position == 0:
                                self.ax1.plot(
                                    self.fh_ax1_factored / qh,
                                    self.fv_ax1_factored / qv,
                                    "--",
                                    color="b",
                                )
                                self.ax1.plot(0, (wbfo - bs) / qv, ".", ms=5, color="b")

        self.handles, self.labels = self.ax1.get_legend_handles_labels()
        self.ax4.legend(
            self.handles,
            self.labels,
            bbox_to_anchor=(0, -0.1),
            loc="upper left",
            prop={"size": self.fontSize},
        )
        self.add_graph_labels()
        self.draw_reacts()
        self.canvas.draw()


if __name__ == "__main__":
    app = App(False)
    frame = VHMEnvelope(parent=None)
    frame.Show()
    app.MainLoop()
