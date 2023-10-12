# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 10:05:52 2020

@author: Wilbert Berendrecht

Class to assign model layer to well screens based on top and bottom
of screens en tops and bots of model layers.
Assumptions:
- Screen length = 1 if one of bkf or okf is missing
- Given the schematization of LHM layers this class uses the top of
  the second layer as bottom of the first layer
- If bkf > top[0] (surface level) and okf < top[0], then screen
  is assigned to layer 1

"""
import os
from hdsrhipy.groundwater.nc import nc
import numpy as np
import configparser


class ModelLayer(object):

    def __init__(self, model='LHM', lhmfolder=None):
        """Class to find model layer in which well screen is located
        """
        config = configparser.ConfigParser()
        config.read('config.ini')

        self.lhmTop = nc(os.path.join(lhmfolder,'TOP.nc'), 'mdl_top')
        self.lhmBot = nc(os.path.join(lhmfolder,'BOT.nc'), 'mdl_bot')


    @staticmethod
    def findTopBot(bkf, okf, bots, tops):
        """Method to find layer index for top and bottom of well screen.

        Parameters
        ----------
        bkf : float
            top of screen (+NAP)
        okf: float
            bot of screen (+NAP)
        x: float
            x-coordinate of screen location
        y: float
            y-coordinate of screen location

        Returns
        -------
        toplaynr: int
            layer in which top of screen is located
        botlaynr: int
            layer in which bot of screen is located
        """
        toplaynr = -1
        botlaynr = -1
        for k in range(len(tops)):
            if toplaynr < 0 and np.isfinite(bots[k]):
                if bkf > bots[k] and np.isfinite(tops[k]):
                    if bkf <= tops[k]:
                        toplaynr = k
                    if okf < tops[k] and okf >= bots[k]:
                        botlaynr = k
            if botlaynr < 0 and toplaynr >=0 and np.isfinite(tops[k]):
                if okf < tops[k] and np.isfinite(bots[k]):
                    if okf >= bots[k]:
                        botlaynr = k

        if toplaynr == -1 and bkf > tops[0] and botlaynr > -1:
            toplaynr = 0

        return toplaynr, botlaynr


    def getLHMLayerRange(self, bkf, okf, x, y):
        """Method to find range of LHM layers in which screen is
           located.

        Parameters
        ----------
        bkf : float
            top of screen (+NAP)
        okf: float
            bot of screen (+NAP)
        x: float
            x-coordinate of screen location
        y: float
            y-coordinate of screen location

        Returns
        -------
        toplaynr: int
            layer in which top of screen is located
        botlaynr: int
            layer in which bot of screen is located
        """

        toplaynr = -1
        botlaynr = -1
        if (np.all(np.isfinite([x, y, bkf, okf])) and
            (x>=self.lhmTop.xmin and
             x<=self.lhmTop.xmax and
             y>=self.lhmTop.ymin and
             y<=self.lhmTop.ymax)
            ):
            tops = self.lhmTop.getValxy(x, y)
            bots = self.lhmBot.getValxy(x, y)
            # due to schematization of layers
            # use top of second layer as bottom of first layer
            bots[0] = tops[1]
            for k in range(len(tops)):
                if toplaynr < 0 and np.isfinite(bots[k]):
                    if bkf >= bots[k]:
                        toplaynr = k
                    elif (k < len(tops) - 1) and (bkf > tops[k+1]):
                        toplaynr = k + 1
            for k in reversed(range(len(tops))):
                if botlaynr < 0 and np.isfinite(tops[k]):
                    if okf <= tops[k]:
                        botlaynr = k
                    elif (k > 0) and (okf < bots[k-1]):
                        botlaynr = k - 1

        return toplaynr, botlaynr


    def getLHMLayer(self, bkf, okf, x, y):
        """Method to find LHM layer number in which
           largest part of screen is located.
        Parameters
        ----------
        bkf : float
            top of screen (+NAP)
        okf: float
            bot of screen (+NAP)
        x: float
            x-coordinate of screen location
        y: float
            y-coordinate of screen location

        Returns
        -------
        lay: int
            layer to which screen is assigned
        """
        lay = np.NaN
        # assume that screen length = 1 if one of bkf or okf is missing
        if np.isnan(okf) and np.isfinite(bkf):
            okf = bkf - 1
        elif np.isnan(bkf) and np.isfinite(okf):
            bkf = okf + 1
        if (np.all(np.isfinite([x,y,bkf,okf])) and
            (bkf - okf <= 10) and
            (x>=self.lhmTop.xmin and
             x<=self.lhmTop.xmax and
             y>=self.lhmTop.ymin and
             y<=self.lhmTop.ymax)
            ):
            tops = self.lhmTop.getValxy(x, y)
            bots = self.lhmBot.getValxy(x, y)
            # due to schematization of layers
            # use top of second layer as bottom of first layer
            bots[0] = tops[1]

            lay = np.NaN
            toplaynr, botlaynr = self.findTopBot(bkf, okf, bots, tops)
            if toplaynr>=0 and botlaynr>=0:
                lay = toplaynr
                # initialize maxLay to store layer with max screen length
                # during search process. Only if toplaynr != botlaynr. If
                # top and bot equal, then whole screen is in single layer
                maxLayThk = -1
                maxLay = -1
                while (toplaynr != botlaynr and toplaynr>=0 and botlaynr>=0):
                    bot1 = bots[toplaynr]
                    len1 = bkf-bot1  # screen length in top layer
                    top2 = tops[botlaynr]
                    len2 = top2-okf # screen length in bot layer
                    if bot1 == top2: # screen is in 2 layers only
                        if len2 > len1:
                            if len2 > maxLayThk:
                                # biggest part of screen is in botlay
                                lay = botlaynr
                            else:
                                # biggest part of screen is in previously
                                # stored maxLay
                                lay = maxLay
                        elif len1 > maxLayThk:
                            # biggest part of screen is in toplay
                            lay = toplaynr
                        else:
                            # biggest part of screen is in previously
                            # stored maxLay
                            lay = maxLay
                        # set toplay = botlay to end while loop
                        toplaynr = botlaynr
                    else:
                        if len2 > len1 and len2 >= (bot1-top2):
                            lay = botlaynr
                            toplaynr=botlaynr
                        elif len1 > len2 and len1 >= (bot1-top2):
                            lay = toplaynr
                            botlaynr=toplaynr
                        elif (bot1-top2) > len1:
                            if len1 > len2 and len1 > maxLayThk:
                                maxLay = toplaynr
                                maxLayThk = len1
                            elif len2 > maxLayThk:
                                maxLay = botlaynr
                                maxLayThk = len2
                            toplaynr, botlaynr = self.findTopBot(bot1, top2,
                                                                bots, tops)
                            if not(toplaynr==-1 and botlaynr==-1):
                                if toplaynr == botlaynr:
                                    lenLay = tops[toplaynr] - bots[toplaynr]
                                    if lenLay >= maxLayThk:
                                        maxLay = toplaynr
                                        maxLayThk = lenLay
                                else:
                                    bkf = bot1
                                    okf = top2
                            else:
                                toplaynr = maxLay
                                botlaynr = maxLay
                            lay = maxLay
            elif botlaynr > 0:
                lay = botlaynr
            elif toplaynr > 0:
                lay = toplaynr

        return lay + 1