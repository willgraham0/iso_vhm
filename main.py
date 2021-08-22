import bisect
import os
import re

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
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
    ALIGN_RIGHT,
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


np.seterr(invalid='ignore', divide='ignore')


class VHMEnvelope(Frame):

    def __init__(self, parent):
        super().__init__(parent, -1, 'ISO V-H-M Envelope ', style=DEFAULT_FRAME_STYLE ^ RESIZE_BORDER)

        # Frame items

        panel = Panel(self)
        panel.SetBackgroundColour('white')

        # Panel items

        self.qv_text = StaticText(panel, -1, 'Qv:')
        self.qv_value = SpinCtrl(panel, -1, min=0, max=100000, initial=10000, size=(75, -1))
        self.qh_text = StaticText(panel, -1, 'Qh:')
        self.qh_value = SpinCtrl(panel, -1, min=0, max=100000, initial=1000, size=(75, -1))
        self.qm_text = StaticText(panel, -1, 'Qm:')
        self.qm_value = SpinCtrl(panel, -1, min=0, max=100000, initial=20000, size=(75, -1))

        self.foundation = RadioBox(panel, -1, label='Foundation', choices=['Sand', 'Clay'])
        self.unit = RadioBox(panel, -1, label='Unit', choices=['Dimensioned', 'Dimensionless'])
        self.suction_box = RadioBox(panel, -1, label='Suction', choices=['Available', 'Unavailable'])
        self.alpha_box = RadioBox(panel, -1, label='Adhesion', choices=[u'\u0251>=0.5', u'\u0251<0.5'])

        self.diameter_text = StaticText(panel, -1, ' Bmax/B: ')
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
        self.a_text = StaticText(panel, -1, ' a ')
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
        self.alpha_text = StaticText(panel, -1, u' \u0251 ')
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

        self.Fv_value_default = self.qv_value.GetValue() / 2
        self.Fh_value_default = 0.0
        self.Fm_value_default = 0.0
        self.Fv_position_text = StaticText(panel, -1, 'Fv:')
        self.Fv_position_value = Slider(
            panel,
            -1,
            value=self.Fv_value_default,
            minValue=0,
            maxValue=self.qv_value.GetValue(),
            size=(150, -1),
            style=SL_LABELS,
        )
        self.Fh_position_text = StaticText(panel, -1, 'Fh:')
        self.Fh_position_value = Slider(
            panel,
            -1,
            value=self.Fh_value_default,
            minValue=0,
            maxValue=self.qh_value.GetValue(),
            size=(150, -1),
            style=SL_LABELS,
        )
        self.Fm_position_text = StaticText(panel, -1, 'Fm:')
        self.Fm_position_value = Slider(
            panel,
            -1,
            value=self.Fm_value_default,
            minValue=0,
            maxValue=self.qm_value.GetValue(),
            size=(150, -1),
            style=SL_LABELS,
        )

        self.draw_button = Button(panel, -1, 'Draw Plots')
        self.clear_button = Button(panel, -1, 'Clear Plots and Reacts')
        self.import_button = Button(panel, -1, 'Import ISO Reactions')
        self.clear_react_button = Button(panel, -1, 'Clear Reactions')
        self.save_button = Button(panel, -1, 'Save Selection')
        self.clear_save_button = Button(panel, -1, 'Clear Selection')
        self.plot_selection_button = Button(panel, -1, 'Plot Selection')

        self.wbfo_value = SpinCtrl(panel, -1, min=0, max=10000, initial=200, size=(150, -1))
        self.wbfo_text = StaticText(panel, -1, 'Wbfo:')
        self.bs_value = SpinCtrl(panel, -1, min=0, max=1000, initial=20, size=(150, -1))
        self.bs_text = StaticText(panel, -1, 'Bs:')
        self.factored_vh = RadioBox(panel, -1, label='Factored V-H', choices=['Off', 'On'])

        # Selection attributes

        self.colours = ['b', 'g', 'r']
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
        self.selection_Fh_ax1 = []
        self.selection_Fv_ax1 = []
        self.selection_Fh_ax2 = []
        self.selection_Fm_ax2 = []
        self.selection_Fv_ax3 = []
        self.selection_Fm_ax3 = []
        self.selection_Fh_ax4 = []
        self.selection_Fv_ax4 = []
        self.selection_Fm_ax4 = []
        self.selection_labels = []

        # Graph attributes

        self.dpi = 100
        self.fontSize = 8

        self.gs = gridspec.GridSpec(3, 2)
        self.fig = plt.figure(facecolor='white', dpi=self.dpi)
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)

        self.ax1 = self.fig.add_subplot(self.gs[0, 0])
        self.ax2 = self.fig.add_subplot(self.gs[1, 0])
        self.ax3 = self.fig.add_subplot(self.gs[2, 0])
        self.ax4 = self.fig.add_subplot(self.gs[0:2, 1], projection='3d')

        # Sizers

        vertical_capacity = BoxSizer(HORIZONTAL)
        vertical_capacity.Add(self.qv_text, flag=ALIGN_CENTRE_VERTICAL)
        vertical_capacity.Add(self.qv_value, flag=ALIGN_CENTRE_VERTICAL)

        horizontal_capacity = BoxSizer(HORIZONTAL)
        horizontal_capacity.Add(self.qh_text, flag=ALIGN_CENTRE_VERTICAL)
        horizontal_capacity.Add(self.qh_value, flag=ALIGN_CENTRE_VERTICAL)

        moment_capacity = BoxSizer(HORIZONTAL)
        moment_capacity.Add(self.qm_text, flag=ALIGN_CENTRE_VERTICAL)
        moment_capacity.Add(self.qm_value, flag=ALIGN_CENTRE_VERTICAL)

        contact_diameter = BoxSizer(HORIZONTAL)
        contact_diameter.Add(self.diameter_text, flag=ALIGN_CENTRE_VERTICAL)
        contact_diameter.Add(self.diameter_value, flag=ALIGN_CENTRE_VERTICAL)

        depth_parameter = BoxSizer(HORIZONTAL)
        depth_parameter.Add(self.a_text, flag=ALIGN_CENTRE_VERTICAL)
        depth_parameter.Add(self.a_value, flag=ALIGN_CENTRE_VERTICAL)

        adhesion = BoxSizer(HORIZONTAL)
        adhesion.Add(self.alpha_text, flag=ALIGN_CENTRE_VERTICAL)
        adhesion.Add(self.alpha_value, flag=ALIGN_CENTRE_VERTICAL)

        capacities = BoxSizer(HORIZONTAL)
        capacities.Add(vertical_capacity)
        capacities.Add(horizontal_capacity)
        capacities.Add(moment_capacity)

        ratios = BoxSizer(HORIZONTAL)
        ratios.Add(contact_diameter)
        ratios.Add(depth_parameter)
        ratios.Add(adhesion)

        fields = BoxSizer(HORIZONTAL)
        fields.Add(self.foundation)
        fields.Add(self.unit)
        fields.Add(self.suction_box)
        fields.Add(self.alpha_box)

        parameters = BoxSizer(HORIZONTAL)
        parameters.Add(capacities, flag=ALIGN_CENTRE_VERTICAL)
        parameters.Add(ratios, flag=ALIGN_CENTRE_VERTICAL)

        sizer_two_sub1 = BoxSizer(HORIZONTAL)
        sizer_two_sub1.Add(self.Fv_position_text, flag=ALIGN_CENTRE_VERTICAL | ALIGN_LEFT)
        sizer_two_sub1.Add(self.Fv_position_value, flag=ALIGN_CENTRE_VERTICAL)

        sizer_two_sub2 = BoxSizer(HORIZONTAL)
        sizer_two_sub2.Add(self.Fh_position_text, flag=ALIGN_CENTRE_VERTICAL | ALIGN_LEFT)
        sizer_two_sub2.Add(self.Fh_position_value, flag=ALIGN_CENTRE_VERTICAL)

        sizer_two_sub3 = BoxSizer(HORIZONTAL)
        sizer_two_sub3.Add(self.Fm_position_text, flag=ALIGN_CENTRE_VERTICAL | ALIGN_LEFT)
        sizer_two_sub3.Add(self.Fm_position_value, flag=ALIGN_CENTRE_VERTICAL)

        sizer_two_sub4 = BoxSizer(HORIZONTAL)
        sizer_two_sub4.Add(self.factored_vh, flag=ALIGN_CENTRE_VERTICAL)  # | ALIGN_RIGHT)

        sizer_two_sub5 = BoxSizer(HORIZONTAL)
        sizer_two_sub5.Add(self.wbfo_text, flag=ALIGN_CENTRE_VERTICAL)  # | ALIGN_RIGHT)
        sizer_two_sub5.Add(self.wbfo_value, flag=ALIGN_CENTRE_VERTICAL)

        sizer_two_sub6 = BoxSizer(HORIZONTAL)
        sizer_two_sub6.Add(self.bs_text, flag=ALIGN_CENTRE_VERTICAL)  # | ALIGN_RIGHT)
        sizer_two_sub6.Add(self.bs_value, flag=ALIGN_CENTRE_VERTICAL)

        sizer_two = BoxSizer(HORIZONTAL)
        sizer_two.Add(sizer_two_sub1)
        sizer_two.Add(sizer_two_sub2)
        sizer_two.Add(sizer_two_sub3)
        sizer_two.Add(sizer_two_sub4, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)
        sizer_two.Add(sizer_two_sub5, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)
        sizer_two.Add(sizer_two_sub6, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)

        actions = BoxSizer(HORIZONTAL)
        actions.Add(self.draw_button, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)
        actions.Add(self.clear_button, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)
        actions.Add(self.import_button, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)
        actions.Add(self.clear_react_button, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)
        actions.Add(self.save_button, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)
        actions.Add(self.clear_save_button, flag=ALIGN_CENTRE_VERTICAL)  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)
        actions.Add(self.plot_selection_button, flag=ALIGN_CENTRE_VERTICAL )  # ALIGN_RIGHT | ALIGN_CENTRE_VERTICAL)

        view = BoxSizer(VERTICAL)
        view.Add(fields)
        view.Add(parameters)
        view.Add(sizer_two, flag=EXPAND)
        view.Add(actions, flag=EXPAND)
        view.Add(self.canvas, 1, flag=LEFT | TOP | GROW)
        view.Add(self.toolbar, flag=EXPAND)

        panel.SetSizer(view)
        view.Fit(self)

        # Events

        self.Bind(EVT_RADIOBOX, self.OnFoundationBox, self.foundation)

        self.Bind(EVT_SPINCTRL, self.OnQvSpin, self.qv_value)
        self.Bind(EVT_SPINCTRL, self.OnQhSpin, self.qh_value)
        self.Bind(EVT_SPINCTRL, self.OnQmSpin, self.qm_value)

        self.Bind(EVT_BUTTON, self.DrawGraphs, self.draw_button)
        self.Bind(EVT_BUTTON, self.ClearGraphs, self.clear_button)
        self.Bind(EVT_BUTTON, self.ImportReacts, self.import_button)
        self.Bind(EVT_BUTTON, self.SaveSelection, self.save_button)
        self.Bind(EVT_BUTTON, self.ClearSelection, self.clear_save_button)
        self.Bind(EVT_BUTTON, self.PlotSelection, self.plot_selection_button)
        self.Bind(EVT_BUTTON, self.ClearReactsWithEvent, self.clear_react_button)

        self.Bind(EVT_RADIOBOX, self.OnFactoredVHBox, self.factored_vh)

        # Initial attribute statuses

        self.a_value.SetValue(0.0)
        self.a_value.Disable()
        self.alpha_value.Disable()
        self.alpha_box.Disable()
        self.suction_box.Disable()
        self.selection_counter = 0
        self.first_plot_flag = 0

        self.wbfo_value.Disable()
        self.bs_value.Disable()

        # Initially executed methods

        self.AddGraphLabels()

        # Constants

        self.gamma_rvh = 1.10

    # Functions

    def get_index(self, list_data, x, flag):
        i = bisect.bisect_left(list_data, x, lo=0, hi=len(list_data))
        if x == list_data[i]:
            return i
        elif x != list_data[i] and flag == 'left':
            return i - 1
        elif x != list_data[i] and flag == 'right':
            return i

    def lookup_table(self, table, value_a, value_alpha):
        a_values = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
        alpha_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        i = self.get_index(a_values, value_a, 'left')
        as_lower = table[:, i]
        j = self.get_index(a_values, value_a, 'right')
        as_upper = table[:, j]

        ii = self.get_index(alpha_values, value_alpha, 'left')
        jj = self.get_index(alpha_values, value_alpha, 'right')

        if i == j and ii == jj:
            value = as_lower[ii]
            return value

        if i == j and ii != jj:
            lower = as_lower[ii]
            upper = as_lower[jj]
            value = (value_alpha - alpha_values[ii]) * (upper - lower) / (alpha_values[jj] - alpha_values[ii]) + lower
            return value

        if i != j and ii == jj:
            lower = as_lower[ii]
            upper = as_upper[ii]
            value = (value_a - a_values[i]) / (a_values[j] - a_values[i]) * (upper - lower) + lower
            return value

        if i != j and ii != jj:
            lower_lower = as_lower[ii]
            lower_upper = as_lower[jj]
            upper_lower = as_upper[ii]
            upper_upper = as_upper[jj]

            alphas_lower = (value_a - a_values[i]) * (as_upper[ii] - as_lower[ii]) / (a_values[j] - a_values[i]) + \
                           as_lower[ii]
            alphas_upper = (value_a - a_values[i]) * (as_upper[jj] - as_lower[jj]) / (a_values[j] - a_values[i]) + \
                           as_lower[jj]
            value = (value_alpha - alpha_values[ii]) * (alphas_upper - alphas_lower) / (
                        alpha_values[jj] - alpha_values[ii]) + alphas_lower
            return value

    def makeFactorisation(self, verticalValues, horizontalValues, wbfo, bs):

        verts = verticalValues - (wbfo - bs)
        vector_lengths = np.sqrt(verts ** 2 + horizontalValues ** 2)
        factored_vector_lengths = (1 / self.gamma_rvh) * vector_lengths

        tan_vector_lengths = verts / vector_lengths
        cos_vector_lengths = horizontalValues / vector_lengths

        factoredVerticals = factored_vector_lengths * tan_vector_lengths + (wbfo - bs)
        factoredHorizontals = factored_vector_lengths * cos_vector_lengths
        return factoredVerticals, factoredHorizontals

    # Methods

    def ClearAxes(self):
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.ax4.clear()

    def PlotSelection(self, event):
        self.ClearAxes()

        for i in range(len(self.selection_Fh_ax1)):
            Qv = self.selection_qv[i]
            Qh = self.selection_qh[i]
            Qm = self.selection_qm[i]

            moments = np.asarray(self.selection_moments[i])
            horizontals = np.asarray(self.selection_horizontals[i])
            verticals = np.asarray(self.selection_verticals[i])

            Fh_ax1 = np.asarray(self.selection_Fh_ax1[i])
            Fv_ax1 = np.asarray(self.selection_Fv_ax1[i])
            Fh_ax2 = np.asarray(self.selection_Fh_ax2[i])
            Fm_ax2 = np.asarray(self.selection_Fm_ax2[i])
            Fv_ax3 = np.asarray(self.selection_Fv_ax3[i])
            Fm_ax3 = np.asarray(self.selection_Fm_ax3[i])
            Fh_ax4 = np.asarray(self.selection_Fh_ax4[i])
            Fv_ax4 = np.asarray(self.selection_Fv_ax4[i])
            Fm_ax4 = np.asarray(self.selection_Fm_ax4[i])
            labels = self.selection_labels[i]

            if self.unit.GetSelection() == 0:
                self.ax1.plot(Fh_ax1, Fv_ax1, color=self.colours[i], label=labels)
                self.ax2.plot(Fh_ax2, Fm_ax2, color=self.colours[i])
                self.ax3.plot(Fv_ax3, Fm_ax3, color=self.colours[i])
                xs, ys = np.meshgrid(Fv_ax4, Fh_ax4)
                self.ax4.plot_surface(xs, ys, Fm_ax4, linewidth=0.1, alpha=0.3, color=self.colours[i])

                self.ax1.scatter(horizontals, verticals, color=self.colours[i], s=5, edgecolor='none')
                self.ax2.scatter(horizontals, moments, color=self.colours[i], s=5, edgecolor='none')
                self.ax3.scatter(verticals, moments, color=self.colours[i], s=5, edgecolor='none')
                self.ax4.scatter(verticals, horizontals, moments, color=self.colours[i], s=3, linewidth='0')

            if self.unit.GetSelection() == 1:
                self.ax1.plot(Fh_ax1 / Qh, Fv_ax1 / Qv, color=self.colours[i], label=labels)
                self.ax2.plot(Fh_ax2 / Qh, Fm_ax2 / Qm, color=self.colours[i])
                self.ax3.plot(Fv_ax3 / Qv, Fm_ax3 / Qm, color=self.colours[i])
                xs, ys = np.meshgrid(Fv_ax4 / Qv, Fh_ax4 / Qh)
                self.ax4.plot_surface(xs, ys, Fm_ax4 / Qm, linewidth=0.1, alpha=0.3, color=self.colours[i])

                self.ax1.scatter(horizontals / Qh, verticals / Qv, color=self.colours[i], s=5, edgecolor='none')
                self.ax2.scatter(horizontals / Qh, moments / Qm, color=self.colours[i], s=5, edgecolor='none')
                self.ax3.scatter(verticals / Qv, moments / Qm, color=self.colours[i], s=5, edgecolor='none')
                self.ax4.scatter(verticals / Qv, horizontals / Qh, moments / Qm, color=self.colours[i], s=3,
                                 linewidth='0')

        self.handles, self.labels = self.ax1.get_legend_handles_labels()
        self.ax4.legend(self.handles, self.labels, bbox_to_anchor=(0, -0.1), loc='upper left',
                        prop={'size': self.fontSize})
        self.AddGraphLabels()
        self.canvas.draw()

    def ClearSelection(self, event):
        self.selection_counter = 0
        self.save_button.Enable()
        self.selection_qv = []
        self.selection_qh = []
        self.selection_qm = []
        self.selection_moments = []
        self.selection_horizontals = []
        self.selection_verticals = []
        self.selection_Fh_ax1 = []
        self.selection_Fv_ax1 = []
        self.selection_Fh_ax2 = []
        self.selection_Fm_ax2 = []
        self.selection_Fv_ax3 = []
        self.selection_Fm_ax3 = []
        self.selection_Fh_ax4 = []
        self.selection_Fv_ax4 = []
        self.selection_Fm_ax4 = []
        self.selection_labels = []

    def ClearReactsWithEvent(self, event):
        self.ClearReacts()

    def SaveSelection(self, event):
        if self.selection_counter <= 2:
            if len(self.moments) == 0:
                self.selection_moments.append([])
                self.selection_horizontals.append([])
                self.selection_verticals.append([])
            else:
                self.selection_moments.append(self.moments)
                self.selection_horizontals.append(self.horizontals)
                self.selection_verticals.append(self.verticals)

            self.selection_qv.append(self.Qv_value.GetValue())
            self.selection_qh.append(self.Qh_value.GetValue())
            self.selection_qm.append(self.Qm_value.GetValue())

            self.selection_Fh_ax1.append(self.Fh_ax1)
            self.selection_Fv_ax1.append(self.Fv_ax1)
            self.selection_Fh_ax2.append(self.Fh_ax2)
            self.selection_Fm_ax2.append(self.Fm_ax2)
            self.selection_Fv_ax3.append(self.Fv_ax3)
            self.selection_Fm_ax3.append(self.Fm_ax3)
            self.selection_Fh_ax4.append(self.Fh_ax4)
            self.selection_Fv_ax4.append(self.Fv_ax4)
            self.selection_Fm_ax4.append(self.Fm_ax4)
            self.selection_labels.append(self.label)

            self.selection_counter = self.selection_counter + 1
            if self.selection_counter > 2:
                self.save_button.Disable()

    def OnFoundationBox(self, event):
        if self.foundation.GetSelection() == 0:
            self.suction_box.Disable()
            self.alpha_box.Disable()
            self.a_value.SetValue(0.00)
            self.alpha_value.SetValue(0.5)
            self.a_value.Disable()
            self.alpha_value.Disable()
            self.diameter_value.SetValue(1.00)
            self.diameter_value.Enable()

        if self.foundation.GetSelection() == 1:
            self.diameter_value.SetValue(1.00)
            self.diameter_value.Disable()
            self.a_value.Enable()
            self.alpha_value.Enable()
            self.alpha_box.Enable()
            self.suction_box.Enable()

    def OnFactoredVHBox(self, event):
        if self.factored_vh.GetSelection() == 0:
            self.wbfo_value.Disable()
            self.bs_value.Disable()

        if self.factored_vh.GetSelection() == 1:
            self.wbfo_value.Enable()
            self.bs_value.Enable()

    def OnQvSpin(self, event):
        self.Fv_position_value.SetMax(self.Qv_value.GetValue())
        self.Fv_position_value.SetValue(self.Qv_value.GetValue() / 2)

    def OnQhSpin(self, event):
        self.Fh_position_value.SetMax(self.Qh_value.GetValue())
        self.Fh_position_value.SetValue(self.Fh_value_default)

    def OnQmSpin(self, event):
        self.Fm_position_value.SetMax(self.Qm_value.GetValue())
        self.Fm_position_value.SetValue(self.Fm_value_default)

    def OnInput(self, event):
        self.DrawGraphs()

    def ImportReacts(self, event):

        dialog = FileDialog(None, 'Choose a .JAR file', defaultDir=os.getcwd(), style=FD_OPEN)
        if dialog.ShowModal() == ID_OK:
            path = dialog.GetPath()
            inFile = open(path, 'r')
            data = inFile.read()
            inFile.close()
            reactions = re.findall(
                r'\d{1,3}\.\d\s*[a-d]\s*\d\s*[\-0]\.\d{4}E[\+\-]\d{2}\s*\-\d*\.\s*\d*\.\s*([\-0]\.\d{4}E[\+\-]\d{2})\s*\-(\d*\.)\s*(\d*\.)\s*\d*\.\s*\d\.\d{2}\s*\d\s*[\-0]\.\d{4}E[\+\-]\d{2}\s*\-\d*\.\s*\d*\.\s*([\-0]\.\d{4}E[\+\-]\d{2})\s*\-(\d*\.)\s*(\d*\.)\s*\d*\.\s*\d\.\d{2}\s*\d\s*[\-0]\.\d{4}E[\+\-]\d{2}\s*\-\d*\.\s*\d*\.\s*([\-0]\.\d{4}E[\+\-]\d{2})\s*\-(\d*\.)\s*(\d*\.)\s*\d*\.\s*\d\.\d{2}',
                data)

            self.moments = np.asarray([[float(i[0]), float(i[3]), float(i[6])] for i in reactions]).flatten()
            self.horizontals = np.asarray([[float(i[1]), float(i[4]), float(i[7])] for i in reactions]).flatten()
            self.verticals = np.asarray([[float(i[2]), float(i[5]), float(i[8])] for i in reactions]).flatten()

    def NormaliseReacts(self):

        self.norm_moments = np.array(self.moments) / float(self.Qm_value.GetValue())
        self.norm_horizontals = np.array(self.horizontals) / float(self.Qh_value.GetValue())
        self.norm_verticals = np.array(self.verticals) / float(self.Qv_value.GetValue())

    def DrawReacts(self):

        self.NormaliseReacts()

        if self.unit.GetSelection() == 0:
            self.ax1.scatter(self.horizontals, self.verticals, s=5, edgecolor='none')
            self.ax2.scatter(self.horizontals, self.moments, s=5, edgecolor='none')
            self.ax3.scatter(self.verticals, self.moments, s=5, edgecolor='none')
            self.ax4.scatter(self.verticals, self.horizontals, self.moments, s=5, linewidth='0')
            self.canvas.draw()

        elif self.unit.GetSelection() == 1:
            self.ax1.scatter(self.norm_horizontals, self.norm_verticals, s=5, edgecolor='none')
            self.ax2.scatter(self.norm_horizontals, self.norm_moments, s=5, edgecolor='none')
            self.ax3.scatter(self.norm_verticals, self.norm_moments, s=5, edgecolor='none')
            self.ax4.scatter(self.norm_verticals, self.norm_horizontals, self.norm_moments, s=5, linewidth='0')
            self.canvas.draw()

        self.AddGraphLabels()

    def AddGraphLabels(self):
        factor = 1.4
        self.ax1.grid(True)
        self.ax2.grid(True)
        self.ax3.grid(True)
        self.ax4.grid(False)

        self.ax1.tick_params(axis='both', which='major', labelsize=self.fontSize)
        self.ax2.tick_params(axis='both', which='major', labelsize=self.fontSize)
        self.ax3.tick_params(axis='both', which='major', labelsize=self.fontSize)
        self.ax4.tick_params(axis='both', which='major', labelsize=self.fontSize)

        if self.unit.GetSelection() == 0:
            self.ax1.set_xlabel(r'$\mathregular{F_H}$ (tonnes)', fontsize=self.fontSize)
            self.ax2.set_xlabel(r'$\mathregular{F_H}$ (tonnes)', fontsize=self.fontSize)
            self.ax3.set_xlabel(r'$\mathregular{F_V}$ (tonnes)', fontsize=self.fontSize)
            self.ax4.set_xlabel(r'$\mathregular{F_V}$ (tonnes)', fontsize=self.fontSize)

            self.ax1.set_ylabel(r'$\mathregular{F_V}$ (tonnes)', fontsize=self.fontSize)
            self.ax2.set_ylabel(r'$\mathregular{F_M}$ (tonnes)', fontsize=self.fontSize)
            self.ax3.set_ylabel(r'$\mathregular{F_M}$ (tonnes)', fontsize=self.fontSize)
            self.ax4.set_ylabel(r'$\mathregular{F_H}$ (tonnes)', fontsize=self.fontSize)

            self.ax4.set_zlabel(r'Moment Reaction, $\mathregular{F_M}$ (tonnes)', fontsize=self.fontSize)

        elif self.unit.GetSelection() == 1:
            self.ax1.set_xlabel(r'$\mathregular{F_H/Q_H}$', fontsize=self.fontSize)
            self.ax2.set_xlabel(r'$\mathregular{F_H/Q_H}$', fontsize=self.fontSize)
            self.ax3.set_xlabel(r'$\mathregular{F_V/Q_V}$', fontsize=self.fontSize)
            self.ax4.set_xlabel(r'$\mathregular{F_V/Q_V}$', fontsize=self.fontSize)

            self.ax1.set_ylabel(r'$\mathregular{F_V/Q_V}$', fontsize=self.fontSize)
            self.ax2.set_ylabel(r'$\mathregular{F_M/Q_M}$', fontsize=self.fontSize)
            self.ax3.set_ylabel(r'$\mathregular{F_M/Q_M}$', fontsize=self.fontSize)
            self.ax4.set_ylabel(r'$\mathregular{F_H/Q_H}$', fontsize=self.fontSize)

            self.ax4.set_zlabel(r'$\mathregular{F_M/Q_M}$', fontsize=self.fontSize)

        self.ax1.set_xlim(left=0)
        self.ax1.set_ylim(bottom=0)
        self.ax2.set_xlim(left=0)
        self.ax2.set_ylim(bottom=0)
        self.ax3.set_xlim(left=0)
        self.ax3.set_ylim(bottom=0)

    def ClearReacts(self):
        self.moments = []
        self.horizontals = []
        self.verticals = []
        self.norm_moments = []
        self.norm_horizontals = []
        self.norm_verticals = []

    def ClearGraphs(self, event):
        self.ClearAxes()
        self.ClearReacts()
        self.AddGraphLabels()
        self.canvas.draw()

    def DrawGraphs(self, event):
        self.ClearAxes()

        FvQvts = np.array([[0.354, 0.334, 0.308, 0.293, 0.276, 0.238, 0.200],
                           [0.387, 0.373, 0.354, 0.343, 0.331, 0.300, 0.265],
                           [0.418, 0.408, 0.396, 0.388, 0.379, 0.357, 0.329],
                           [0.447, 0.441, 0.433, 0.428, 0.423, 0.409, 0.390],
                           [0.474, 0.471, 0.468, 0.465, 0.463, 0.457, 0.448],
                           [0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500]])

        m_alphas = np.array([[1.172, 1.200, 1.239, 1.264, 1.295, 1.378, 1.500],
                             [0.902, 0.917, 0.937, 0.950, 0.965, 1.006, 1.067],
                             [0.653, 0.661, 0.670, 0.676, 0.683, 0.701, 0.729],
                             [0.422, 0.425, 0.429, 0.431, 0.434, 0.440, 0.450],
                             [0.205, 0.206, 0.207, 0.207, 0.208, 0.209, 0.211],
                             [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]])

        granularity = 300
        Qv = float(self.qv_value.GetValue())
        Qh = float(self.qh_value.GetValue())
        Qm = float(self.qm_value.GetValue())
        Bm_B = float(self.diameter_value.GetValue())
        a = float(self.a_value.GetValue())
        alpha = float(self.alpha_value.GetValue())
        alpha_flag = self.alpha_box.GetSelection()
        Fv_position = float(self.Fv_position_value.GetValue())
        Fh_position = float(self.Fh_position_value.GetValue())
        Fm_position = float(self.Fm_position_value.GetValue())
        wbfo = float(self.wbfo_value.GetValue())
        bs = float(self.bs_value.GetValue())

        if self.foundation.GetSelection() == 0:
            soil_flag = 'Sand'

        elif self.foundation.GetSelection() == 1:
            soil_flag = 'Clay'
            suction_flag = self.suction_box.GetSelection()

        if self.unit.GetSelection() == 0:
            graph_flag = 'Dimensioned'

        elif self.unit.GetSelection() == 1:
            graph_flag = 'Dimensionless'

        startv = 0.0
        starth = 0.0
        FvQvt = self.lookup_table(FvQvts, a, alpha)
        m_alpha = self.lookup_table(m_alphas, a, alpha)

        if soil_flag == 'Sand':
            Fv = np.linspace(startv, Qv, granularity)

            self.Fv_ax1 = np.linspace(startv, Qv, granularity)
            Qmps_ax1 = Qm * (Bm_B) ** 3 * np.ones(granularity)
            Qmpv_ax1 = 0.15 * ((Qm * 0.12) / (0.075 * Qh)) * self.Fv_ax1
            Fh_ax1_orig = np.nan_to_num(Qh * np.sqrt(
                4 * a * (self.Fv_ax1 / Qv) * (1 - (self.Fv_ax1 / Qv)) + 16 * (1 - a) * (self.Fv_ax1 / Qv) ** 2 * (
                            1 - (self.Fv_ax1 / Qv)) ** 2 - (Fm_position / Qm) ** 2))
            Fh_ax1_mod = np.nan_to_num(Qh * np.sqrt(
                4 * a * (self.Fv_ax1 / Qv) * (1 - (self.Fv_ax1 / Qv)) + 16 * (1 - a) * (self.Fv_ax1 / Qv) ** 2 * (
                            1 - (self.Fv_ax1 / Qv)) ** 2 - (Fm_position / np.minimum(Qmps_ax1, Qmpv_ax1)) ** 2))
            self.Fh_ax1 = np.maximum(Fh_ax1_orig, Fh_ax1_mod)

            if Fm_position == 0:
                self.Fv_ax1_factored, self.Fh_ax1_factored = self.makeFactorisation(self.Fv_ax1, self.Fh_ax1, wbfo, bs)

            Fv_ax2 = np.linspace(startv, Qv, granularity)
            self.Fh_ax2 = np.linspace(starth, Qh, granularity)
            Qmps_ax2 = Qm * (Bm_B) ** 3 * np.ones(granularity)
            Qmpv_ax2 = 0.15 * ((Qm * 0.12) / (0.075 * Qh)) * Fv_position
            Fm_ax2_orig = np.nan_to_num(Qm * np.sqrt(
                4 * a * (Fv_position / Qv) * (1 - (Fv_position / Qv)) + 16 * (1 - a) * (Fv_position / Qv) ** 2 * (
                            1 - (Fv_position / Qv)) ** 2 - (self.Fh_ax2 / Qh) ** 2))
            Fm_ax2_mod = np.nan_to_num(np.minimum(Qmps_ax2, Qmpv_ax2) * np.sqrt(
                4 * a * (Fv_position / Qv) * (1 - (Fv_position / Qv)) + 16 * (1 - a) * (Fv_position / Qv) ** 2 * (
                            1 - (Fv_position / Qv)) ** 2 - (self.Fh_ax2 / Qh) ** 2))
            self.Fm_ax2 = np.maximum(Fm_ax2_orig, Fm_ax2_mod)

            self.Fv_ax3 = np.linspace(startv, Qv, granularity)
            Qmps_ax3 = Qm * (Bm_B) ** 3 * np.ones(granularity)
            Qmpv_ax3 = 0.15 * ((Qm * 0.12) / (0.075 * Qh)) * self.Fv_ax3
            Fm_ax3_orig = Qm * np.sqrt(
                4 * a * (self.Fv_ax3 / Qv) * (1 - (self.Fv_ax3 / Qv)) + 16 * (1 - a) * (self.Fv_ax3 / Qv) ** 2 * (
                            1 - (self.Fv_ax3 / Qv)) ** 2 - (Fh_position / Qh) ** 2)
            Fm_ax3_mod = np.minimum(Qmps_ax3, Qmpv_ax3) * np.sqrt(
                4 * a * (self.Fv_ax3 / Qv) * (1 - (self.Fv_ax3 / Qv)) + 16 * (1 - a) * (self.Fv_ax3 / Qv) ** 2 * (
                            1 - (self.Fv_ax3 / Qv)) ** 2 - (Fh_position / Qh) ** 2)
            self.Fm_ax3 = np.maximum(Fm_ax3_orig, Fm_ax3_mod)

            self.Fv_ax4 = np.linspace(startv, Qv, granularity)
            self.Fh_ax4 = np.linspace(starth, Qh, granularity)
            Qmps_ax4 = Qm * (Bm_B) ** 3 * np.ones(granularity)
            Qmpv_ax4 = 0.15 * ((Qm * 0.12) / (0.075 * Qh)) * self.Fv_ax4
            Fm_ax4_before = []

            for i in range(granularity):
                Fm_inter_orig = Qm * np.sqrt(
                    4 * a * (self.Fv_ax4 / Qv) * (1 - self.Fv_ax4 / Qv) + 16 * (1 - a) * (self.Fv_ax4 / Qv) ** 2 * (
                                1 - self.Fv_ax4 / Qv) ** 2 - (self.Fh_ax4[i] / Qh) ** 2)
                Fm_inter_mod = np.minimum(Qmps_ax4, Qmpv_ax4) * np.sqrt(
                    4 * a * (self.Fv_ax4 / Qv) * (1 - self.Fv_ax4 / Qv) + 16 * (1 - a) * (self.Fv_ax4 / Qv) ** 2 * (
                                1 - self.Fv_ax4 / Qv) ** 2 - (self.Fh_ax4[i] / Qh) ** 2)
                Fm_inter = np.maximum(Fm_inter_orig, Fm_inter_mod)
                Fm_ax4_before.append(Fm_inter)

            self.Fm_ax4 = np.nan_to_num(Fm_ax4_before)
            xs, ys = np.meshgrid(self.Fv_ax4, self.Fh_ax4)

            self.label = '$\mathregular{Q_V}$=%dt; $\mathregular{Q_H}$=%dt; $\mathregular{Q_M}$=%dt; $\mathregular{B_{max}}$/B= %.2f' % (
            Qv, Qh, Qm, Bm_B)
            if graph_flag == 'Dimensioned':
                self.ax1.plot(self.Fh_ax1, self.Fv_ax1, label=self.label)
                self.ax2.plot(self.Fh_ax2, self.Fm_ax2)
                self.ax3.plot(self.Fv_ax3, self.Fm_ax3)
                self.ax4.plot_surface(xs, ys, self.Fm_ax4, linewidth=0.1, alpha=0.3)

                if self.factored_vh.GetSelection() == 1:
                    if Fm_position == 0:
                        self.ax1.plot(self.Fh_ax1_factored, self.Fv_ax1_factored, '--', color='b')
                        self.ax1.plot(0, wbfo - bs, '.', ms=5, color='b')
            ##                    if any(self.Fv_ax1_factored):
            ##                        zeds = np.array(len(self.Fh_ax1_factored))
            ##                        self.ax4.plot(self.Fv_ax1_factored, self.Fh_ax1_factored, zeds, '--', color='b')

            elif graph_flag == 'Dimensionless':
                self.ax1.plot(self.Fh_ax1 / Qh, self.Fv_ax1 / Qv, label=self.label)
                self.ax2.plot(self.Fh_ax2 / Qh, self.Fm_ax2 / Qm)
                self.ax3.plot(self.Fv_ax3 / Qv, self.Fm_ax3 / Qm)
                self.ax4.plot_surface(xs / Qv, ys / Qh, self.Fm_ax4 / Qm, linewidth=0.1, alpha=0.3)

                if self.factored_vh.GetSelection() == 1:
                    if Fm_position == 0:
                        self.ax1.plot(self.Fh_ax1_factored / Qh, self.Fv_ax1_factored / Qv, '--', color='b')
                        self.ax1.plot(0, (wbfo - bs) / Qv, '.', ms=5, color='b')

        if soil_flag == 'Clay':
            if alpha_flag == 1:
                Fv = np.linspace(startv, Qv, granularity)

                self.Fv_ax1 = np.linspace(startv, Qv, granularity)
                self.Fh_ax1 = Qh * np.sqrt(
                    4 * a * (self.Fv_ax1 / Qv) * (1 - (self.Fv_ax1 / Qv)) + 16 * (1 - a) * (self.Fv_ax1 / Qv) ** 2 * (
                                1 - (self.Fv_ax1 / Qv)) ** 2 - (Fm_position / Qm) ** 2)

                if Fm_position == 0:
                    self.Fv_ax1_factored, self.Fh_ax1_factored = self.makeFactorisation(self.Fv_ax1, self.Fh_ax1, wbfo,
                                                                                        bs)

                self.Fh_ax2 = np.linspace(starth, Qh, granularity)
                self.Fm_ax2 = Qm * np.sqrt(
                    4 * a * (Fv_position / Qv) * (1 - (Fv_position / Qv)) + 16 * (1 - a) * (Fv_position / Qv) ** 2 * (
                                1 - (Fv_position / Qv)) ** 2 - (self.Fh_ax2 / Qh) ** 2)

                self.Fv_ax3 = np.linspace(startv, Qv, granularity)
                self.Fm_ax3 = Qm * np.sqrt(
                    4 * a * (self.Fv_ax3 / Qv) * (1 - (self.Fv_ax3 / Qv)) + 16 * (1 - a) * (self.Fv_ax3 / Qv) ** 2 * (
                                1 - (self.Fv_ax3 / Qv)) ** 2 - (Fh_position / Qh) ** 2)

                self.Fv_ax4 = np.linspace(startv, Qv, granularity)
                self.Fh_ax4 = np.linspace(starth, Qh, granularity)
                Fm_ax4_before = []

                for i in range(granularity):
                    Fm_inter = Qm * np.sqrt(
                        4 * a * (self.Fv_ax4 / Qv) * (1 - self.Fv_ax4 / Qv) + 16 * (1 - a) * (self.Fv_ax4 / Qv) ** 2 * (
                                    1 - self.Fv_ax4 / Qv) ** 2 - (self.Fh_ax4[i] / Qh) ** 2)
                    Fm_ax4_before.append(Fm_inter)
                self.Fm_ax4 = np.nan_to_num(Fm_ax4_before)
                xs, ys = np.meshgrid(self.Fv_ax4, self.Fh_ax4)

                self.label = '$\mathregular{Q_V}$=%dt; $\mathregular{Q_H}$=%dt; $\mathregular{Q_M}$=%dt; $\mathregular{\\alpha}$<0.50; a= %.2f; Suction N/A' % (
                Qv, Qh, Qm, a)
                if graph_flag == 'Dimensioned':
                    self.ax1.plot(self.Fh_ax1, self.Fv_ax1, label=self.label)
                    self.ax2.plot(self.Fh_ax2, self.Fm_ax2)
                    self.ax3.plot(self.Fv_ax3, self.Fm_ax3)
                    self.ax4.plot_surface(xs, ys, self.Fm_ax4, linewidth=0.1, alpha=0.3)

                    if self.factored_vh.GetSelection() == 1:
                        if Fm_position == 0:
                            self.ax1.plot(self.Fh_ax1_factored, self.Fv_ax1_factored, '--', color='b')
                            self.ax1.plot(0, wbfo - bs, '.', ms=5, color='b')

                elif graph_flag == 'Dimensionless':
                    self.ax1.plot(self.Fh_ax1 / Qh, self.Fv_ax1 / Qv, label=self.label)
                    self.ax2.plot(self.Fh_ax2 / Qh, self.Fm_ax2 / Qm)
                    self.ax3.plot(self.Fv_ax3 / Qv, self.Fm_ax3 / Qm)
                    self.ax4.plot_surface(xs / Qv, ys / Qh, self.Fm_ax4 / Qm, linewidth=0.1, alpha=0.3)

                if self.factored_vh.GetSelection() == 1:
                    if Fm_position == 0:
                        self.ax1.plot(self.Fh_ax1_factored / Qh, self.Fv_ax1_factored / Qv, '--', color='b')
                        self.ax1.plot(0, (wbfo - bs) / Qv, '.', ms=5, color='b')

            if alpha_flag == 0:
                if suction_flag == 0:
                    Fv = np.linspace(Qv * FvQvt, Qv, granularity)

                    Fv_ax1_part1 = Qv * np.linspace(0.0, FvQvt, granularity)
                    Fv_ax1_part2 = np.linspace(Qv * FvQvt, Qv, granularity)
                    self.Fv_ax1 = np.append(Fv_ax1_part1, Fv_ax1_part2)
                    Fh_ax1_part1 = Qh * np.sqrt((alpha + m_alpha * (Fv_ax1_part1 / Qv)) ** 2 - (Fm_position / Qm) ** 2)
                    Fh_ax1_part2 = Qh * np.sqrt(
                        4 * a * (Fv_ax1_part2 / Qv) * (1 - (Fv_ax1_part2 / Qv)) + 16 * (1 - a) * (
                                    Fv_ax1_part2 / Qv) ** 2 * (1 - (Fv_ax1_part2 / Qv)) ** 2 - (Fm_position / Qm) ** 2)
                    self.Fh_ax1 = np.append(Fh_ax1_part1, Fh_ax1_part2)

                    if Fm_position == 0:
                        self.Fv_ax1_factored, self.Fh_ax1_factored = self.makeFactorisation(self.Fv_ax1, self.Fh_ax1,
                                                                                            wbfo, bs)

                    Fh_part1_lim = Qh * np.sqrt((alpha + m_alpha * (FvQvt)) ** 2 - (0 / Qm) ** 2)
                    Fh_ax2_part1 = np.linspace(0.0, Fh_part1_lim, granularity)
                    Fh_ax2_part2 = np.linspace(0.0, Qh, granularity)
                    Fm_ax2_part1 = Qm * np.sqrt((alpha + m_alpha * (Fv_position / Qv)) ** 2 - (Fh_ax2_part1 / Qh) ** 2)
                    Fm_ax2_part2 = Qm * np.sqrt(4 * a * (Fv_position / Qv) * (1 - (Fv_position / Qv)) + 16 * (1 - a) * (
                                Fv_position / Qv) ** 2 * (1 - (Fv_position / Qv)) ** 2 - (Fh_ax2_part2 / Qh) ** 2)

                    if Fv_position <= Qv * FvQvt:
                        self.Fh_ax2 = Fh_ax2_part1
                        self.Fm_ax2 = Fm_ax2_part1
                    elif Fv_position > Qv * FvQvt:
                        self.Fh_ax2 = Fh_ax2_part2
                        self.Fm_ax2 = Fm_ax2_part2

                    Fv_ax3_part1 = Qv * np.linspace(0.0, FvQvt, granularity)
                    Fv_ax3_part2 = np.linspace(Qv * FvQvt, Qv, granularity)
                    self.Fv_ax3 = np.append(Fv_ax3_part1, Fv_ax3_part2)
                    Fm_ax3_part1 = Qm * np.sqrt((alpha + m_alpha * (Fv_ax3_part1 / Qv)) ** 2 - (Fh_position / Qh) ** 2)
                    Fm_ax3_part2 = Qm * np.sqrt(
                        4 * a * (Fv_ax3_part2 / Qv) * (1 - (Fv_ax3_part2 / Qv)) + 16 * (1 - a) * (
                                    Fv_ax3_part2 / Qv) ** 2 * (1 - (Fv_ax3_part2 / Qv)) ** 2 - (Fh_position / Qh) ** 2)
                    self.Fm_ax3 = np.append(Fm_ax3_part1, Fm_ax3_part2)

                    Fv_ax4_part1 = np.linspace(0.0, Qv * FvQvt, granularity / 2)
                    self.Fh_ax4 = np.linspace(0.0, Qh, granularity)
                    Fm_ax4_part1 = []

                    for i in range(granularity):
                        Fm_ax4_part1_inter = Qm * np.sqrt(
                            (alpha + m_alpha * (Fv_ax4_part1 / Qv)) ** 2 - (self.Fh_ax4[i] / Qh) ** 2)
                        Fm_ax4_part1.append(Fm_ax4_part1_inter)
                    Fm_ax4_part1_nan = np.nan_to_num(Fm_ax4_part1)

                    Fv_ax4_part2 = np.linspace(Qv * FvQvt, Qv, granularity / 2)
                    Fm_ax4_part2 = []

                    for i in range(granularity):
                        Fm_ax4_part2_inter = Qm * np.sqrt(
                            4 * a * (Fv_ax4_part2 / Qv) * (1 - Fv_ax4_part2 / Qv) + 16 * (1 - a) * (
                                        Fv_ax4_part2 / Qv) ** 2 * (1 - Fv_ax4_part2 / Qv) ** 2 - (
                                        self.Fh_ax4[i] / Qh) ** 2)
                        Fm_ax4_part2.append(Fm_ax4_part2_inter)
                    Fm_ax4_part2_nan = np.nan_to_num(Fm_ax4_part2)

                    self.Fv_ax4 = np.append(Fv_ax4_part1, Fv_ax4_part2)
                    self.Fm_ax4 = np.append(Fm_ax4_part1_nan, Fm_ax4_part2_nan, axis=1)
                    Xs, Ys = np.meshgrid(self.Fv_ax4, self.Fh_ax4)

                    self.label = '$\mathregular{Q_V}$=%dt; $\mathregular{Q_H}$=%dt; $\mathregular{Q_M}$=%dt; $\mathregular{\\alpha}$= %.2f; a= %.2f; Suction' % (
                    Qv, Qh, Qm, alpha, a)
                    if graph_flag == 'Dimensioned':
                        self.ax1.plot(self.Fh_ax1, self.Fv_ax1, label=self.label)
                        self.ax2.plot(self.Fh_ax2, self.Fm_ax2)
                        self.ax3.plot(self.Fv_ax3, self.Fm_ax3)
                        self.ax4.plot_surface(Xs, Ys, self.Fm_ax4, linewidth=0.1, alpha=0.3)

                        if self.factored_vh.GetSelection() == 1:
                            if Fm_position == 0:
                                self.ax1.plot(self.Fh_ax1_factored, self.Fv_ax1_factored, '--', color='b')
                                self.ax1.plot(0, wbfo - bs, '.', ms=5, color='b')

                    elif graph_flag == 'Dimensionless':
                        self.ax1.plot(self.Fh_ax1 / Qh, self.Fv_ax1 / Qv, label=self.label)
                        self.ax2.plot(self.Fh_ax2 / Qh, self.Fm_ax2 / Qm)
                        self.ax3.plot(self.Fv_ax3 / Qv, self.Fm_ax3 / Qm)
                        self.ax4.plot_surface(Xs / Qv, Ys / Qh, self.Fm_ax4 / Qm, linewidth=0.1, alpha=0.3)

                        if self.factored_vh.GetSelection() == 1:
                            if Fm_position == 0:
                                self.ax1.plot(self.Fh_ax1_factored / Qh, self.Fv_ax1_factored / Qv, '--', color='b')
                                self.ax1.plot(0, (wbfo - bs) / Qv, '.', ms=5, color='b')

                if suction_flag == 1:

                    Fv = np.linspace(Qv * FvQvt, Qv, granularity)

                    Fv_ax1_part1 = np.linspace(1, Qv * FvQvt, granularity)
                    Fv_ax1_part2 = np.linspace(Qv * FvQvt, Qv, granularity)
                    self.Fv_ax1 = np.append(Fv_ax1_part1, Fv_ax1_part2)
                    Fh_ax1_part1 = Qh * (alpha + m_alpha * (Fv_ax1_part1 / Qv)) * np.sqrt(1 - (Fm_position / (
                                Qm * np.sqrt(
                            16 * (1 - a) * (Fv_ax1_part1 / Qv) ** 2 * (1 - (Fv_ax1_part1 / Qv)) ** 2 + 4 * a * (
                                        Fv_ax1_part1 / Qv) * (1 - (Fv_ax1_part1 / Qv))))) ** 2)
                    Fh_ax1_part2 = Qh * np.sqrt(
                        4 * a * (Fv_ax1_part2 / Qv) * (1 - (Fv_ax1_part2 / Qv)) + 16 * (1 - a) * (
                                    Fv_ax1_part2 / Qv) ** 2 * (1 - (Fv_ax1_part2 / Qv)) ** 2 - (Fm_position / Qm) ** 2)
                    self.Fh_ax1 = np.append(Fh_ax1_part1, Fh_ax1_part2)

                    if Fm_position == 0:
                        self.Fv_ax1_factored, self.Fh_ax1_factored = self.makeFactorisation(self.Fv_ax1, self.Fh_ax1,
                                                                                            wbfo, bs)

                    Fh_part1_lim = Qh * np.sqrt((alpha + m_alpha * (FvQvt)) ** 2 - (0 / Qm) ** 2)
                    Fh_ax2_part1 = np.linspace(0.0, Fh_part1_lim, granularity)
                    Fh_ax2_part2 = np.linspace(0.0, Qh, granularity)
                    Fm_ax2_part1 = (Qm * np.sqrt(
                        16 * (1 - a) * (Fv_position / Qv) ** 2 * (1 - (Fv_position / Qv)) ** 2 + 4 * a * (
                                    Fv_position / Qv) * (1 - (Fv_position / Qv)))) * np.sqrt(
                        1 - (Fh_ax2_part1 / ((alpha + m_alpha * (Fv_position / Qv)) * Qh)) ** 2)
                    Fm_ax2_part2 = Qm * np.sqrt(4 * a * (Fv_position / Qv) * (1 - (Fv_position / Qv)) + 16 * (1 - a) * (
                                Fv_position / Qv) ** 2 * (1 - (Fv_position / Qv)) ** 2 - (Fh_ax2_part2 / Qh) ** 2)

                    if Fv_position <= Qv * FvQvt:
                        self.Fh_ax2 = Fh_ax2_part1
                        self.Fm_ax2 = Fm_ax2_part1
                    elif Fv_position > Qv * FvQvt:
                        self.Fh_ax2 = Fh_ax2_part2
                        self.Fm_ax2 = Fm_ax2_part2

                    Fv_ax3_part1 = Qv * np.linspace(0.0, FvQvt, granularity)
                    Fv_ax3_part2 = np.linspace(Qv * FvQvt, Qv, granularity)
                    self.Fv_ax3 = np.append(Fv_ax3_part1, Fv_ax3_part2)
                    Fm_ax3_part1 = (Qm * np.sqrt(
                        16 * (1 - a) * (Fv_ax3_part1 / Qv) ** 2 * (1 - (Fv_ax3_part1 / Qv)) ** 2 + 4 * a * (
                                    Fv_ax3_part1 / Qv) * (1 - (Fv_ax3_part1 / Qv)))) * np.sqrt(
                        1 - (Fh_position / ((alpha + m_alpha * (Fv_ax3_part1 / Qv)) * Qh)) ** 2)
                    Fm_ax3_part2 = Qm * np.sqrt(
                        4 * a * (Fv_ax3_part2 / Qv) * (1 - (Fv_ax3_part2 / Qv)) + 16 * (1 - a) * (
                                    Fv_ax3_part2 / Qv) ** 2 * (1 - (Fv_ax3_part2 / Qv)) ** 2 - (Fh_position / Qh) ** 2)
                    self.Fm_ax3 = np.append(Fm_ax3_part1, Fm_ax3_part2)

                    Fv_ax4_part1 = np.linspace(0.0, Qv * FvQvt, granularity / 2)
                    self.Fh_ax4 = np.linspace(0.0, Qh, granularity)
                    Fm_ax4_part1 = []

                    for i in range(granularity):
                        Fm_ax4_part1_inter = (Qm * np.sqrt(
                            16 * (1 - a) * (Fv_ax4_part1 / Qv) ** 2 * (1 - (Fv_ax4_part1 / Qv)) ** 2 + 4 * a * (
                                        Fv_ax4_part1 / Qv) * (1 - (Fv_ax4_part1 / Qv)))) * np.sqrt(
                            1 - (self.Fh_ax4[i] / ((alpha + m_alpha * (Fv_ax4_part1 / Qv)) * Qh)) ** 2)
                        Fm_ax4_part1.append(Fm_ax4_part1_inter)
                    Fm_ax4_part1_nan = np.nan_to_num(Fm_ax4_part1)

                    Fv_ax4_part2 = np.linspace(Qv * FvQvt, Qv, granularity / 2)
                    Fm_ax4_part2 = []

                    for i in range(granularity):
                        Fm_ax4_part2_inter = Qm * np.sqrt(
                            4 * a * (Fv_ax4_part2 / Qv) * (1 - Fv_ax4_part2 / Qv) + 16 * (1 - a) * (
                                        Fv_ax4_part2 / Qv) ** 2 * (1 - Fv_ax4_part2 / Qv) ** 2 - (
                                        self.Fh_ax4[i] / Qh) ** 2)
                        Fm_ax4_part2.append(Fm_ax4_part2_inter)
                    Fm_ax4_part2_nan = np.nan_to_num(Fm_ax4_part2)

                    self.Fv_ax4 = np.append(Fv_ax4_part1, Fv_ax4_part2)
                    self.Fm_ax4 = np.append(Fm_ax4_part1_nan, Fm_ax4_part2_nan, axis=1)
                    Xs, Ys = np.meshgrid(self.Fv_ax4, self.Fh_ax4)

                    self.label = '$\mathregular{Q_V}$=%dt; $\mathregular{Q_H}$=%dt; $\mathregular{Q_M}$=%dt; $\mathregular{\\alpha}$= %.2f; a= %.2f; No Suction' % (
                    Qv, Qh, Qm, alpha, a)
                    if graph_flag == 'Dimensioned':
                        self.ax1.plot(self.Fh_ax1, self.Fv_ax1, label=self.label)
                        self.ax2.plot(self.Fh_ax2, self.Fm_ax2)
                        self.ax3.plot(self.Fv_ax3, self.Fm_ax3)
                        self.ax4.plot_surface(Xs, Ys, self.Fm_ax4, linewidth=0.1, alpha=0.3)

                        if self.factored_vh.GetSelection() == 1:
                            if Fm_position == 0:
                                self.ax1.plot(self.Fh_ax1_factored, self.Fv_ax1_factored, '--', color='b')
                                self.ax1.plot(0, wbfo - bs, '.', ms=5, color='b')

                    elif graph_flag == 'Dimensionless':
                        self.ax1.plot(self.Fh_ax1 / Qh, self.Fv_ax1 / Qv, label=self.label)
                        self.ax2.plot(self.Fh_ax2 / Qh, self.Fm_ax2 / Qm)
                        self.ax3.plot(self.Fv_ax3 / Qv, self.Fm_ax3 / Qm)
                        self.ax4.plot_surface(Xs / Qv, Ys / Qh, self.Fm_ax4 / Qm, linewidth=0.1, alpha=0.3)

                        if self.factored_vh.GetSelection() == 1:
                            if Fm_position == 0:
                                self.ax1.plot(self.Fh_ax1_factored / Qh, self.Fv_ax1_factored / Qv, '--', color='b')
                                self.ax1.plot(0, (wbfo - bs) / Qv, '.', ms=5, color='b')

        self.handles, self.labels = self.ax1.get_legend_handles_labels()
        self.ax4.legend(self.handles, self.labels, bbox_to_anchor=(0, -0.1), loc='upper left',
                        prop={'size': self.fontSize})
        self.AddGraphLabels()
        self.DrawReacts()
        self.canvas.draw()


if __name__ == '__main__':
    app = App(False)
    frame = VHMEnvelope(parent=None)
    frame.Show()
    app.MainLoop()
