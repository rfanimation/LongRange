#!/usr/bin/python

import csv
import math


def findAlphaOmega(tableFilename):
    """ This looks up the range value at the top and bottom of the
        ballistic table. It is used wherever the difference in range
        values (hi to low vs. low to high) must be considered.
         Inputs:  tableFilename = the table to be read
         Outputs: none
         returns: tableTop and table Bottom
    """
    import os
    import sys
    scriptPath = os.path.dirname(sys.argv[0])
    tableFilename = os.path.join(scriptPath, 'ingalls.csv')
    outFilename = os.path.join(scriptPath, 'tableOut.csv')
    reader = csv.reader(open(tableFilename, "rb"))
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        tableTop = float(row[0])
        break
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        tableBottom = float(row[0])
    return tableTop, tableBottom

def calcDrop(v, V, tof):
    """ calculate bullet drop from V/Vo and time of flight
        No interpolation, use closest value
         Inputs:
           v = velocity in fps at range
           V = muzzle velocity in fps
           tof = time of flight in seconds
           drop.csv = drop table
         Outputs: none
         returns: drop in inches
        Note 1: drop table from Hatcher's Notebook
        Note to programmers: The drop table is only used here
          and there is no planned alternative table so the name
          is hard coded rather than passed purposefully.
    """
    tgt = v/V       # target match value
    last1 = 1.0
    import os
    import sys
    scriptPath = os.path.dirname(sys.argv[0])
    dropFilename = os.path.join(scriptPath, 'drop.csv')
    reader = csv.reader(open(dropFilename, "rb"))
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        # adjust so closest vale V/Vo will be used
        adjusted = float(row[0]) - 0.005
        if adjusted < tgt:
            f = float(row[1])
            retVal = (tof * tof) * f
            return retVal
    # not found error exit
    errVal = (tof * tof) * 137.0
    #print "V/Vo is too small, drop too large to calculate."
    return 0.0

def interpolate(r, inhi, inlo, outhi, outlo):
    """ Performs a simple linear interpolation
         Inputs:
           r = reference value (typically range)
           inhi = hi value for r value
           inlo = lo value for r value
           outhi = hi value to be scaled
           outlo = lo value to be scaled
         Outputs: none
         returns: scaled value from out column
    """
    # calculate scale of reference value
    tmp01 = inhi-inlo
    if 0.0 == tmp01:
        tmp01 = 0.000000000001
    m = (-1*(r-inhi)/tmp01)
    # scale the output value
    retVal = outlo + ((outhi - outlo) * m)
    return retVal

def lookupST4range(u, tableFilename, tableTop):
    """ Fetches space time data for range in the form of
         [ u, S(u), delta, T(u), delta ]
         Inputs:
           u = range to be looked up
           tableFilename = name of table file
           tableTop = used to determine table type
         Outputs: none
         returns: space time data
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
        Note 2: If interpolation is required only S(u) and
        T(u) are interpolated as deltas are not required for
        any function fetching this data.
        Note 3: If the look-up goes beyond the table limits
          an obviously bogus -1 is returned to flag the fact
          that the value could not be calculated.
    """
    import os
    import sys
    scriptPath = os.path.dirname(sys.argv[0])
    tableFilename = os.path.join(scriptPath, 'ingalls.csv')
    outFilename = os.path.join(scriptPath, 'tableOut.csv')

    reader = csv.reader(open(tableFilename, "rb"))
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        # 1st look for an exact match
        if u == float(row[0]):
            retVal = row
            # convert text to numeric
            retVal[4] = float(row[4])
            retVal[3] = float(row[3])
            retVal[2] = float(row[2])
            retVal[1] = float(row[1])
            retVal[0] = float(row[0])
            return retVal
    # No exact match, must interpolate
    reader = csv.reader(open(tableFilename, "rb"))
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        lastRow = row
        lastU = float(row[0])
        break
    reader = csv.reader(open(tableFilename, "rb"))
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        if tableTop < tableBottom:
            # British or other slow to fast table
            if float(row[0]) > u:
                retVal = row
                # interpolate Su
                retVal[1] = interpolate(u, float(row[0]), float(lastRow[0]), float(row[1]), float(lastRow[1]))
                # interpolate T(u)
                retVal[3] = interpolate(u, float(row[0]), float(lastRow[0]), float(row[3]), float(lastRow[3]))
                # put the range in the output
                retVal[0] = u
                retVal[2] = float(row[2]) # delta not interpolated but converted for table checking
                retVal[4] = float(row[4]) # delta not interpolated but converted for table checking
                # return range and interpolated S(u) and T(u). Deltas unchanged
                return retVal
        else:
            # Ingalls or other fast to slow table
            if float(row[0]) < u:
                retVal = row
                # interpolate Su
                retVal[1] = interpolate(u, float(lastRow[0]), float(row[0]), float(row[1]), float(lastRow[1]))
                # interpolate T(u)
                retVal[3] = interpolate(u, float(lastRow[0]), float(row[0]), float(row[3]), float(lastRow[3]))
                # put the range in the output
                retVal[0] = u
                retVal[2] = float(row[2]) # delta not interpolated but converted for table checking
                retVal[4] = float(row[4]) # delta not interpolated but converted for table checking
                # return range and interpolated S(u) and T(u). Deltas unchanged
                return retVal
        lastU = float(row[0])
        lastRow = row
    print "lookupST4range(",u,",",tableFilename,",",tableTop,") failed"
    retVal = [-1,-1,-1,-1,-1]
    return retVal

def lookuprange4ST(SV, tableFilename, tableTop):
    """ Fetches range for Sv from space time data table
         Inputs:
           Sv = space value from table
           tableFilename = name of table file
           tableTop = used to determine table type
         Outputs: none
         returns: range in feet
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
    """
    import os
    import sys
    scriptPath = os.path.dirname(sys.argv[0])
    tableFilename = os.path.join(scriptPath, 'ingalls.csv')
    outFilename = os.path.join(scriptPath, 'tableOut.csv')
    reader = csv.reader(open(tableFilename, "rb"))
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        # 1st look for an exact match
        if SV == float(row[1]):
            return float(row[0])
    # No exact match, must interpolate
    reader = csv.reader(open(tableFilename, "rb"))
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        lastRow = row
        lastSV = float(row[1])
        break
    reader = csv.reader(open(tableFilename, "rb"))
    for row in reader:
        # skip text header
        try:
            xyz = float(row[0])
        except ValueError:
            continue
        if float(row[1]) > SV:
            # Interpolate range
            if tableTop < tableBottom:
                # British or other slow to fast table
                retVal = interpolate(SV, float(row[1]), float(lastRow[1]), float(row[0]), float(lastRow[0]))
                return retVal
            else:
                # Ingalls or other fast to slow table
                retVal = interpolate(SV, float(row[1]), float(lastRow[1]), float(lastRow[0]), float(row[0]))
                return retVal
        lastU = float(row[0])
        lastSV = float(row[1])
        lastRow = row
    print "lookuprange4ST(",SV,",",tableFilename,",",tableTop,") failed"
    #return float(row[0])

def calcCfromChron(V, v, X, tableFilename, tableTop):
    """ Calculates ballistic coefficient (C) from V, v and X
        This facilitates calculation of C when it is not known
        using chronograph measurements.
         Inputs:
           V = muzzle velocity
           v = velocity at distance
           X = distance in feet
           tableFilename = name of table file
           tableTop = used to determine table type
         Outputs: none
         returns: ballistic coefficient
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
    """
    data1 = lookupST4range(V, tableFilename, tableTop)
    data2 = lookupST4range(v, tableFilename, tableTop)
    SV = float(data1[1])
    Sv = float(data2[1])
    C = X /  (Sv - SV)
    return C

def calcFormFactor(V, v, X, w, d, C):
    """ Calculates form factor
        This facilitates calculation of i when it is not known
        using chronograph measurements from C=id<squared>
         Inputs:
           V = initial velocity in fps
           v = velocity at distance
           X = distance in feet
           w = weight in grains
           d = diameter in inches
           tableFilename = name of table file (see note 4)
           BC = previously calculated (see note 4)
         Outputs: none
         returns: form factor
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
        Note 2: w is defined as weight in pounds but we enter the
          weight in grains to this function which converts to pounds.
        Note 3: d is generally taken as bore diameter + depth of one
          rifling groove (effective diameter).
        Note 4: currently, this only is called if calcCfromChron
          has already been called. The input 'C' it's normalised
          result (form factor at standard environment)
    """
    # calculate w/d<squared> given weight in grains
    temp1 = w / (d * d * 7000.0)
    # calculate form factor
    i = temp1 / C
    return i

def calcRemainingVelocity(V, C, X, tableFilename, tableTop):
    """ Calculates remaining velocity (v) at distance X
        Formula: S(v)=S(V)+(X/C)
         Inputs:
           V = initial velocity in fps
           C = ballistic coefficient
           X = distance in feet
           tableFilename = name of table file
           tableTop = used to determine table type
         Outputs: none
         returns: velocity at distance X
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
    """
    data1 = lookupST4range(V, tableFilename, tableTop)
    SV = data1[1]
    if tableTop < tableBottom:
        # British or other slow to fast table
        Sv = SV - (X / C)
    else:
        # Ingalls or other fast to slow table
        Sv = SV + (X / C)
    vel = lookuprange4ST(Sv, tableFilename, tableTop)
    return vel

def calcTimeOfFlight(V, C, X, tableFilename, tableTop):
    """ Calculates time of flight (T) to distance X from
        Formula: T=C * {T(v)-(T(v)}
         Inputs:
           V = initial velocity in fps
           C = ballistic coefficient
           X = distance to time
           tableFilename = name of table file
           tableTop = used to determine table type
         Outputs: none
         returns: time of flight in seconds
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
    """
    data1 = lookupST4range(V, tableFilename, tableTop)
    TV = float(data1[3])
    v = calcRemainingVelocity(V, C, X, tableFilename, tableTop)
    data2 = lookupST4range(v, tableFilename, tableTop)
    Tv = float(data2[3])
    T = (Tv - TV) * C
    if tableTop < tableBottom:
        # British or other slow to fast table
        T = T * -1
    return T

def calcMaxOrdinate(V, C, X, tableFilename):
    """ Calculates maximum ordinate for range X
        Formula:
           H(feet)=(2T)<squared>=4T<squared> H(inches)=48T<squared>
         Inputs:
           V = initial velocity in fps
           C = ballistic coefficient
           X = distance that firearm is sighted in for
           tableFilename = name of table file
         Outputs: none
         returns: maximum ordinate height
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
        Note 2: many people think of this as the midrange trajectory
          but the maximum ordinate is actually a bit farther than
          the midrange point
    """
    T = calcTimeOfFlight(V, C, X, tableFilename, tableTop)
    H = (T * T) * 48.0
    return H

def calcMaxOrdinateRange(X):
    """ Estimates the maximum ordinate distance from range.
         Inputs: X = distance that firearm is sighted in for
         Outputs: none
         returns: maximum ordinate distance
        Note 1: This is an attempt to compensate for the 'negligible'
          difference to mid-range with a magic fudge factor.
          It turns out that this did not improve things. I am leaving
          the function but output value is zero.
    """
    yds = X / 3.0
    temp1 = yds / 2.0
#   temp2 = temp1 * 0.055
    temp2 = temp1 * 0.0
    mor = (temp1 + temp2) * 3.0
    #print yds, temp1, temp2, mor, (mor/yds) * 100
    return mor

def calcAngleOfDeparture(V, C, X, tableFilename):
    """ Calculates angle of departure for sight-in range.
         Inputs:
           V = initial velocity in fps
           C = ballistic coefficient
           X = range the firearm is sighted in for
           tableFilename = name of table file
         Outputs:
         returns:
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
        Note 2: Hatcher's Notebook has a lookup table for this
          which is derived from Ingalls Tables (without extensions)
          which would limit the function to a subset of possible
          ranges. This function uses a method of my own creation.
    """
    H = calcMaxOrdinate(V, C, X, tableFilename)
    mor = calcMaxOrdinateRange(X)
    # prevent divide by zero error
    if 0.0 == mor:
        mor = 1.0e-9
    #print "fixed divide by zero error in calcAngleOfDeparture(", V, ",", C, ",", X, ",", tableFilename, ")"
    v = calcRemainingVelocity(V, C, mor, tableFilename, tableTop)
    T = calcTimeOfFlight(V, C, mor, tableFilename, tableTop)
    drop = calcDrop(v, V, T)
    if 0.0 == drop:
        # return error value, can't calculate drop
        return 0.0
    lod = H + drop
    tanRad = lod / (mor * 12.0)
    moa = (tanRad * 180 / math.pi) * 60.0
    return moa

def calcTemperatureCorrectionFactor(temperature, sTemp):
    """ Calculates a correction factor to ballistic coefficient
        for temperature deviation from standard
        Formula:
          (shooting site temp + 459.4) / (standard temp + 459.4)
         Inputs:  temperature in degrees Fahrenheit
         Outputs: none
         returns: scaling factor
        Note 1: The source of this formula is the Hornady website.
    """
    tFactor = (temperature + 459.4)  / (sTemp + 459.4)
    return tFactor

def calcPressureCorrectionFactor(pressure, sPres):
    """ Calculates a correction factor to ballistic coefficient
        for pressure deviation from standard
        Formula:
          standard pressure / shooting site pressure
         Inputs:  pressure inHg
         Outputs: none
         returns: scaling factor
        Note 1: The source of this formula is the Hornady website.
    """
    pFactor = sPres / pressure
    return pFactor

def calcPressureTempCorrection(BC, pFactor, tFactor):
    """ Returns corrected ballistic coefficient from
        calculated correction factors
         Inputs:
           BC = ballistic coefficient
           pFactor = scaling factor for pressure
           tFactor = scaling factor for temperature
         Outputs: none
         returns: compensated ballistic coefficient
        Note 1: if only one correction factor has been calculated
        enter 1.0 (unity) for the one not calculated.
        Note 2: The source of this formula is the Hornady website.
    """
    C = BC * pFactor * tFactor
    return C

def calcWindDisplacement(W, WA, V, C, X, tableFilename):
    """ Calculates wind displacement in inches
        Formula: D = W*(T-Tv)*F
         Where T = TOF in air, sec. & Tv = TOF in a vacuum, sec
         F = A factor relating to the wind angle.
         Inputs:
           W = Crosswind velocity, fps
           WA = Wind angle in degrees
           V = velocity in fps
           C = ballistic coefficient
           X = range in feet
           tableFilename = name of table file
         Outputs: none
         returns: Bullet displacement, ft.
        Note 1: table file default is Ingalls with extension
          but can be another table with the same format.
        Note 2: formula is for displacement in feet so the
          return value is scaled in the calculation for D
    """
    temp1 = math.radians(WA)
    F = math.sin(temp1)
    temp2 = W * 1.4667
    Wm = temp2
    T = calcTimeOfFlight(V, C, X, tableFilename, tableTop)
    Tv = X / V
    D = Wm * (T - Tv) * F * 12.0
    return D

def calcLead(tof, X, sot):
    """ Calculates lead for moving target in mils
        Formula: lead = tof / (17.6 / inch per mil)
         Inputs:
           tof = time of flight in seconds
           sot = speed of target in MPH
           X = range in yards
         Outputs: none
         returns: lead in mils.
    """
    ipm = (X * 3 * 12) * (1.0 / 1000)
    lead = tof / (17.6 / ipm)
    return lead

def calcEnergy(v, w):
    """ Calculates impact energy in ft/lbs from weight grains and velocity
        Formula: E=1/2mv<squared>=(w*v<squared>/2g
         Inputs:
           V = velocity in fps
           w = weight in grains
         Outputs: none
         returns: energy in foot pounds
         Note 1: formula is for weight in pounds so this function
           scales for use of grains weight.
    """
    E = (w * v * v) / (2.0 * 32.16 * 7000.0)
    return E

def calcImpact(V, v, X, zero, C, sightHeight):
    """ calculates impact with respect to sight picture
         Inputs:
           V = muzzle velocity in fps
           v = velocity at distance in fps
           X = distance in feet
           x = range the firearm is sighted in for
           C = ballistic coefficient
         Outputs: none
         returns: impact above or below sight picture in inches
    """
    if 0.0 == zero:
        zero = 1.0e-9
        #print "fixed divide by zero error in calcImpact(", V, ",", v, ",", X, ",", zero, ",", C, ",", sightHeight, ") "
    # calculate drop
    tof = calcTimeOfFlight(V, C, X, tableFilename, tableTop)
    drop = calcDrop(v, V, tof)
    if 0.0 == drop:
        # return error value, can't calculate drop
        return 0.0
    # calculate height of line of departure above range
    aod = calcAngleOfDeparture(V, C, zero, tableFilename)
    aodDeg = aod / 60.0
    tanRad = math.radians(aodDeg)
    opposite = tanRad * (X * 12.0)
    # calculate distance of line of sight
    aosTan = sightHeight / zero
    sightO = aosTan * (zero - X)
    # calculate the impact distance
    impact = opposite - drop - sightO
    return impact

def addRow(data, newRow):
    """ Build the table with data for a new range
         Inputs:
           data = the existing table
           newrow = the new row of data
         Outputs: updated data
         returns: none
        Note to programmers: Python is normally object oriented. Since
          I intended for do-it-yourselfers to be able to experiment with
          this, I chose not to make it object oriented. A side effect of
          that decision is that data, rather than pointers like in C are
          passed. I pass data to functions, and return it as if that were
          not the case to make my intentions clear. For this function to
          work properly I had to make an unbound copy of newRow using the
          lines lines that precede the append. Don't change that.
    """
    # make an unbounded copy of the list
    tempRow = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    tempRow += newRow
    # discard original nulls
    for letter in "123456789abcdef":
        del tempRow[0]
    data.append(tempRow)

def calculateHit(step, MV, wt, C, temp, pres, zero, sightHeight, sTemp, sPres, sot, range):
    """ Perform calculations for a new row of data for output table.
         Inputs:
         Outputs: Row of data displayed.
         returns: a new row of calculated data
    """

    import csv
    import math

    # convert to feet
    range = range * 3.0
    #newRow[2] = calcEnergy(MV, wt)
    #newRow[7] = sightHeight * -1

    # do temperature and pressure compensation
    if sTemp <> temp:
        tFactor = calcTemperatureCorrectionFactor(temp, sTemp)
    else:
        tFactor = 1.0
    if sPres <> pres:
        pFactor = calcPressureCorrectionFactor(pres, sPres)
    else:
        pFactor = 1.0
    BC = calcPressureTempCorrection(C, pFactor, tFactor)
    C = BC

    # convert to feet
    zeroRange = zero * 3.0
    v = calcRemainingVelocity(MV, C, range, tableFilename, tableTop)
    tof = calcTimeOfFlight(MV, C, range, tableFilename, tableTop)
    mo = calcMaxOrdinate(MV, C, range, tableFilename)
    moRng = calcMaxOrdinateRange(range)
    moMil = inchToMil(mo, range)
    moMin = milToMoa(moMil)
    saod = calcAngleOfDeparture(MV, C, range, tableFilename)
    aod = calcAngleOfDeparture(MV, C, zeroRange, tableFilename)
    drop = calcDrop(v, MV, tof)
    impact = calcImpact(MV, v, range, zeroRange, C, sightHeight)
    iMil = inchToMil(impact, range)
    iMin = milToMoa(iMil)
    drift = calcWindDisplacement(W, WA, MV, C, range, tableFilename)
    dMil = inchToMil(drift, range)
    dMin = milToMoa(dMil)
    E = calcEnergy(v, wt)

    return drop, drift, tof

def inchToMoa(inch, feet):
    """ This function us used to convert inches to moa.
         Inputs:
           inch = opposite side in inches
           feet = adjacent side in feet
         Outputs: none
         returns: angle in minutes
    """
    mil = (inch / (feet * 12)) * 1000
    return milToMoa(mil)

def inchToMil(inch, feet):
    """ This function us used to convert inches to mils
         Inputs:
           inch = opposite side in inches
           feet = adjacent side in feet
         Outputs: none
         returns: angle in milliradians
    """
    # prevent divide by zero error
    if 0.0 == feet:
        feet = 1.0e-9
        #print "fixed divide by zero error in inchToMil(", inch, ",", feet, ")"
    if inch < 0.0:
        inch = inch * -1
        retVal = (inch / (feet * 12)) * 1000
        retVal = retVal * -1
    else:
        retVal = (inch / (feet * 12)) * 1000
    return retVal

def moaToInch(moa, feet):
    """ This function us used to convert moa to inches
         Inputs:
           moa =  angle in minutes
           feet = adjacent side in feet
         Outputs: none
         returns: adjacent side in inches
    """
    mil = moaToMil(moa)
    return milToInch(mil, feet)

def milToInch(mil, feet):
    """ This function us used to convert mils to inches
         Inputs:
           mil =  angle in milliradians
           feet = adjacent side in feet
         Outputs: none
         returns: adjacent side in inches
    """
    return (feet * 12) * (mil / 1000)

def moaToMil(moa):
    """ This function us used
         Inputs: moa = angle in minutes
         Outputs: none
         returns: angle in milliradians
    """
    return moa * 0.29088821

def milToMoa(mil):
    """ This function us used to convert mils to moa
         Inputs: angle in milliradians
         Outputs: none
         returns: angle in minutes
    """
    return mil * 3.4377468

#####################
#  Main starts here #
#####################

#
# initial defaults
#

showAll = 0     # 1=echo input (debug)
helpOn = 0      # 1=on

tableFilename = "ingalls.csv"
outFilename = "tableOut.csv"
tableTop = 0
tableBottom = 0

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
tFactor = calcTemperatureCorrectionFactor(temp, sTemp)
pFactor = calcPressureCorrectionFactor(pres, sPres)
C = calcPressureTempCorrection(BC, pFactor, tFactor)

####
#### START
####

tableTop, tableBottom = findAlphaOmega(tableFilename)

initRng = step * -1.0

range = 1000

drop, drift, tof = calculateHit(step, MV, wt, C, temp, pres, zero, sightHeight, sTemp, sPres, sot, range)
print str(drop)
print str(drift)
print str(tof)


