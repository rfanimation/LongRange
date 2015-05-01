# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
# Name:        LongRange
# Purpose:
#
# Author:      Rob Field
#
# Created:     28/04/2015
# Copyright:   (c) Rob Field 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import time

from kivy.app import App
from kivy.lang import Builder
from kivy.uix.image import Image
from kivy.core.window import Window
from kivy.core.audio import SoundLoader
from kivy.uix.screenmanager import ScreenManager, Screen

class Manager(ScreenManager):
    pass

class TitleScreen(Screen):
    def loadRangeScreen(self):
        self.manager.current = 'range'
    def loadSettingsScreen(self):
        self.manager.current = 'settings'

class RangeScreen(Screen):
    def createTargetPressed(self):
        window_x = Window.width
        window_y = Window.height

        scaleMultiplier = 1.0
        self.steelTarget.size_hint = (scaleMultiplier, scaleMultiplier)
        self.steelTarget.pos_hint={'center_x': .5 , 'center_y': .5}

    def fireButtonPressed(self):
        import csv
        import math
        import ballistics

        #####################
        #  Main starts here #
        #####################

        #
        # initial defaults
        #

        showAll = 0     # 1=echo input (debug)
        helpOn = 0      # 1=on

        # Ballistic coefficient compensation data
        sMetroAlti = 0.0        # sea level
        sMetroTemp = 59.0       # degrees F
        sMetroBaro = 29.5275    # inches Hg
        sMetrohumi = 78.0       # percent
        sICAOAlti = 0.0         # sea level
        sICAOTemp = 59.0        # degrees F
        sICAOBaro = 29.9213     # inches Hg
        sICAOhumi = 0.0         # percent

        sTemp = sMetroTemp      # standard temp
        sPres = sMetroBaro      # standard pressure

        temp = sMetroTemp       # default standard Metro for no compensation
        pres = sMetroBaro       # default standard Metro for no compensation

        tFactor = 1.0   # gets calculated
        pFactor = 1.0   # gets calculated

        #
        # wind data for the chart (drift)
        #
        W = 10.0        # cross wind mph
        WA = 90.0       # cross wind direction

        #
        # bullet data
        #
        BC = 0.4        # published ballistic coefficient
        C = 0.0         # compensated ballistic coefficient (calculated
        wt = 150.0      # projectile weight in grains
        MV = 2700.0     # muzzle velocity in fps
        diam = 0.308    # diameter of projectile

        #
        # for impact calculations
        #
        sightHeight = 1.5    # default scope height

        #
        # calculating distances for output
        #
        step = 50.0     # step interval
        maxRng = 1000.0  # max range
        zero = 250.0    # zero range


        maxErr = 3.0    # max distance from point of aim

        menuItem = 1
        doMore = 1

        sot = 0.0 # Speed of moving target
        calcPBR = 0 # Don't calculate point blank range

        #
        # Compensated (local) ballistic coefficient
        #
        tFactor = ballistics.calcTemperatureCorrectionFactor(temp, sTemp)
        pFactor = ballistics.calcPressureCorrectionFactor(pres, sPres)
        C = ballistics.calcPressureTempCorrection(BC, pFactor, tFactor)

        ####
        #### START
        ####

        initRng = step * -1.0

        range = 100

        drop, drift, tof = ballistics.calculateHit(step, MV, wt, C, temp, pres, zero, sightHeight, sTemp, sPres, sot, range)
        shotSound = SoundLoader.load('shot.wav')
        hitSound = SoundLoader.load('hit.wav')
        shotSound.play()
        time.sleep(tof)
        hitSound.play()
        window_x = Window.width
        window_y = Window.height
        reticle_x = (self.reticleScatter.pos[0] / window_x)
        reticle_y = (self.reticleScatter.pos[1] / window_y)
        self.hitPoint.pos_hint = {'x': reticle_x, 'y': reticle_y}

    def zoomInButtonPressed(self):
        global zoomLevel
        if zoomLevel <= 2:
            zoomLevel = zoomLevel * 2
            if ffpActive == False:
                self.rangeScatter.scale *= 2
            if ffpActive == True:
                self.ffpScatter.scale *= 2

    def zoomOutButtonPressed(self):
        global zoomLevel
        if zoomLevel >= 2:
            zoomLevel = zoomLevel / 2
            if ffpActive == False:
                self.rangeScatter.scale /= 2
            if ffpActive == True:
                self.ffpScatter.scale /= 2

class SettingsScreen(Screen):
    def reloadTitleScreen(self):
        self.manager.current = 'title'

    def toggleButtonPressed(self):
        global ffpActive
        if self.ffpToggle.state == 'down':
            ffpActive = True
        elif self.sfpToggle.state == 'down':
            ffpActive = False
        else:
            self.ffpToggle.state = 'down'
            ffpActive = True


root_widget = Builder.load_string('''

#:kivy 1.0
Manager:
    TitleScreen:
    RangeScreen:
    SettingsScreen:

<TitleScreen>:

    titleLabel: titleLabelID
    rangeScreenButton: rangeScreenButtonID
    settingsScreenButton: settingsScreenButtonID

    name: 'title'

    BoxLayout:
        orientation: 'vertical'
        Label:
            id: titleLabelID
            text: "Long Range"
            text_size: (self.size[0] * 0.8, self.size[1])
            halign: 'center'
            valign: 'middle'
            font_size: self.height / 3

        Button:
            id: rangeScreenButtonID
            text: 'Go to the Range'
            on_release: root.loadRangeScreen()

        Button:
            id: settingsScreenButtonID
            text: 'Settings'
            on_release: root.loadSettingsScreen()

<RangeScreen>:

    ffpScatter: ffpScatterID
    rangeScatter: rangeScatterID
    rangeBackground: rangeBackgroundID
    steelTarget: steelTargetID
    reticleScatter: reticleScatterID
    reticleImage: reticleImageID
    hitPoint: hitPointID
    backToMenu: backToMenuID
    fireButton: fireButtonID
    zoomInButton: zoomInButtonID
    zoomOutButton: zoomOutButtonID
    createTarget: createTargetID

    name: 'range'

    FloatLayout:
        id: box1ID
        orientation: 'vertical'

        ScatterLayout:
            id: ffpScatterID
            do_rotation: False
            do_scale: False
            do_translation: False
            auto_bring_to_front: False
            scale: 1

            ScatterLayout:
                id: rangeScatterID
                do_rotation: False
                do_scale: False
                do_translation: False
                auto_bring_to_front: False
                scale: 1

                Image:
                    id: rangeBackgroundID
                    source: 'rifleRange.jpg'
                    size_hint: (2, 2)
                    pos_hint: {'center_x': 0.465, 'center_y': 0.4}

                Image:
                    id: steelTargetID
                    source: 'ipsc_steel_10.png'


            ScatterLayout:
                id: reticleScatterID
                do_rotation: False
                do_scale: False
                auto_bring_to_front: False
                scale: 1

                Image:
                    id: reticleImageID
                    source: 'milDotReticle.png'

            Image:
                id: hitPointID
                source: 'hitPoint.png'
                pos_hint: {'x': -10000000, 'y': 10000000}

        Button:
            id: backToMenuID
            text: 'Menu'
            text_size: (self.size[0] * 0.8, self.size[1])
            halign: 'center'
            valign: 'middle'
            font_size: self.height / 3
            on_release: app.root.current = 'title'
            size_hint: (0.1, 0.1)

        Button:
            id: fireButtonID
            background_color: (1,0,0,1)
            text: 'FIRE'
            text_size: self.size
            halign: 'center'
            valign: 'middle'
            font_size: self.height / 2
            on_release: root.fireButtonPressed()
            size_hint: (0.1, 0.1)
            pos_hint: {'center_x': 0.5, 'y': 0.0}

        Button:
            id: zoomInButtonID
            text: 'Zoom In'
            text_size: (self.size[0] * 0.8, self.size[1])
            halign: 'center'
            valign: 'middle'
            font_size: self.height / 3
            on_release: root.zoomInButtonPressed()
            size_hint: (0.1, 0.1)
            pos_hint: {'x': 0.0, 'y': 0.9}

        Button:
            id: zoomOutButtonID
            text: 'Zoom Out'
            text_size: (self.size[0] * 0.8, self.size[1])
            halign: 'center'
            valign: 'middle'
            font_size: self.height / 3
            on_release: root.zoomOutButtonPressed()
            size_hint: (0.1, 0.1)
            pos_hint: {'x': 0.0, 'y': 0.8}

        Button:
            id: createTargetID
            text: 'Create Target'
            text_size: (self.size[0] * 0.8, self.size[1])
            halign: 'center'
            valign: 'middle'
            font_size: self.height / 3
            on_release: root.createTargetPressed()
            size_hint: (0.1, 0.1)
            pos_hint: {'x': 0.0, 'y': 0.6}


<SettingsScreen>:

    settingsTitleLabel: settingsTitleLabelID
    titleButton: titleButtonID
    ffpToggle: ffpToggleID
    sfpToggle: sfpToggleID

    name: 'settings'

    BoxLayout:
        orientation: 'vertical'

        Label:
            id: settingsTitleLabelID
            text: "Settings"
            text_size: (self.size[0] * 0.8, self.size[1])
            halign: 'center'
            valign: 'middle'
            font_size: self.height / 3

        Button:
            id: titleButtonID
            text: 'Return to Title Screen'
            on_release: root.reloadTitleScreen()

        BoxLayout:
            orientation: 'horizontal'

            ToggleButton:
                id: ffpToggleID
                text: 'First Focal Plane'
                group: 'focalPlane'
                on_release: root.toggleButtonPressed()

            ToggleButton:
                id: sfpToggleID
                text: 'Second Focal Plane'
                group: 'focalPlane'
                state: 'down'
                on_release: root.toggleButtonPressed()

''')


class LongRange(App):
    def build(self):
        global ffpActive
        ffpActive = False
        global zoomLevel
        zoomLevel = 1
        return root_widget

LongRange().run()
