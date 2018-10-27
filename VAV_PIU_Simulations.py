from random import uniform, random, choice
import math
from datetime import datetime, timedelta

# create the log directory if there is no log directory
import os
simulationLogDirectory = os.path.join(".", "BuildingLogs")
if not os.path.exists(simulationLogDirectory):
    os.makedirs(simulationLogDirectory)

def deleteSimulationLogs(logStartsWith=None):
    if logStartsWith is None:
        # delete out all of the files
        for the_file in os.listdir(simulationLogDirectory):
            file_path = os.path.join(simulationLogDirectory, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                # elif os.path.isdir(file_path): shutil.rmtree(file_path) # performs a recursive delete on any folders
            except Exception as e:
                print(e)
    else:
        # delete out all of the files
        prefixLength = len(logStartsWith)
        for the_file in os.listdir(simulationLogDirectory):
            filePrefix = the_file[:prefixLength]
            if filePrefix == logStartsWith:
                file_path = os.path.join(simulationLogDirectory, the_file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                    # elif os.path.isdir(file_path): shutil.rmtree(file_path) # performs a recursive delete on any folders
                except Exception as e:
                    print(e)

def selectRandomFromRange(startEl, endEl):
    if startEl == endEl:
        return startEl
    return uniform(startEl, endEl)

def selectRandomElement(elList, elLikelihoods=None):
    if elLikelihoods == None:
        return choice(elList)
    else:
        listLen = len(elList)
        likelyhoodLen = len(elLikelihoods)
        assert listLen == likelyhoodLen, "Element list and likelyhood tuples must be the same length."
        assert listLen > 0, "Element list and likelyhood tuples must have length greater than 0."
        randNum = random()
        cumulativeAmount = 0
        for i in range(listLen):
            cumulativeAmount += elLikelihoods[i]
            if randNum < cumulativeAmount:
                return elList[i]
        return elList[-1]  # for if the probabilities don't add up to 1 (or more than 1)

def selectRandomFromRangeList(rangeList, rangeLikelihoods):
    selectedRange = selectRandomElement(rangeList, rangeLikelihoods)
    assert len(selectedRange) == 2, "Range tuple must have length two."
    return selectRandomFromRange(selectedRange[0], selectedRange[1])

def degreesToRadians(degreeAngle):
    return math.pi / 180. * degreeAngle

def radiansToDegrees(radianAngle):
    return 180. / math.pi * radianAngle

def roundToNearest(number, roundToNumber=1.):
    return round(number / roundToNumber, 0) * roundToNumber

def roundUpToNearest(number, roundToNumber=1.):
    return math.ceil(number / roundToNumber) * roundToNumber

def returnSign(number):
    if number < 0.:
        return -1
    else:
        return 1

def solveForExteriorTempOfWall(solarInsolation, wallEmissivity, convectionCoefficient, wallResistance, innerTemperature, outsideTemperature, guessTemperature=None, threshold=0.001):
    # solarInsolation is the heat provided by the sun, in BTU / (hr * ft^2)
    # wallEmmissivity in the radiation emmission constant of the wall, which is between 0 and 1
    # convectionCoefficient is the coefficient of convection, measured in BTU / (hr * ft^2 * degF), due to wind on the outside side of the wall
    # wallResistance is the effective resistance of the wall, in (hr * degF) / BTU
    # innerTemperature is the temperature on the inside of the room (on the other side of the wall) in degF
    # outsideTemperature is the ambient temperature outside of the room in degF
    # guessTemperature is the starting guess in degF
    # threshold is a stopping criteria - stop after the change in guesses is lower than this
    if solarInsolation <= 0:
        totalResistance = 1. / convectionCoefficient + wallResistance
        # high value for convectionCoefficient means low convection resistance
        # if wall resistance >> convection resistance, should be close to outside temperature
        # if wall resistance << convection resistance, should be close to inside temperature
        return innerTemperature + (wallResistance / totalResistance) * (outsideTemperature - innerTemperature)
    else:
        if guessTemperature == None:
            guessTemperature = outsideTemperature
        # convert temperatures to Rankine for usage
        innerTempRankine = innerTemperature + 460
        outerTempRankine = outsideTemperature + 460
        guessTempRankine = guessTemperature + 460

        # set up constants
        stefanBoltzmannConstant = 0.1714e-8  # in BTU / (hr * ft^2 * degF^2)
        quarticCoefficient = stefanBoltzmannConstant * wallEmissivity
        linearCoefficient = convectionCoefficient + 1. / wallResistance
        constantCoefficient = - stefanBoltzmannConstant * wallEmissivity * outerTempRankine ** 4 - convectionCoefficient * outerTempRankine - 1. / wallResistance * outerTempRankine - solarInsolation

        # set up the initial guesses
        point1Temp = guessTempRankine
        point2Temp = guessTempRankine + 20
        point1Value = quarticCoefficient * point1Temp ** 4 + linearCoefficient * point1Temp + constantCoefficient
        point2Value = quarticCoefficient * point2Temp ** 4 + linearCoefficient * point2Temp + constantCoefficient
        currDiff = point2Temp - point1Temp

        # run the loop
        while abs(currDiff) > threshold:
            guessTemp = point1Temp + (point2Temp - point1Temp) * - point1Value / (point2Value - point1Value)
            guessValue = quarticCoefficient * guessTemp ** 4 + linearCoefficient * guessTemp + constantCoefficient
            if abs(point1Value) > abs(point2Value):
                # replace point 1
                point1Value = guessValue
                point1Temp = guessTemp
            else:
                # replace point 2
                point2Value = guessValue
                point2Temp = guessTemp
            currDiff = point2Temp - point1Temp

        return point1Temp - 460



class ConstructionMaterialValues:
    # class which contains information about construction materials, like insulation values and reflectance values for
    # the roof, exterior wall, and interior wall
    def __init__(self):
        self.generateRoofParameters()
        self.generateExteriorWallParameters()
        self.generateExteriorWindowParameters()
        self.generateInteriorWallParameters()
        self.generateFloorParameters()

    def generateRoofParameters(self):
        self.roofRValue = round(selectRandomFromRange(12, 17), 1)
        self.roofEmittance = 0.8  # just made up this number - maybe should randomly generate it?

    def generateExteriorWallParameters(self):
        self.outerWallRValue = round(selectRandomFromRange(12, 17), 1)
        self.outerWallEmittance = 0.8  # just made up this number - maybe should randomly generate it?

    def generateExteriorWindowParameters(self):
        self.windowRValue = round(selectRandomFromRange(2, 4), 1)
        self.windowTransmittance = selectRandomElement((0.25, 0.4, 0.6), (0.4, 0.4, 0.2))

    def generateInteriorWallParameters(self):
        self.innerWallRValue = round(selectRandomFromRange(3, 5), 1)

    def generateFloorParameters(self):
        self.floorRValue = round(selectRandomFromRange(7, 12), 1)

class AngleHandler:
    # class to handle 3-dimensional angles, such as sun angles, wind directions, and wall normal faces
    def __init__(self, angleFromHorizontal=0., azimuthAngle=0., vectorLength=1.):
        self.setAngles(angleFromHorizontal, azimuthAngle)
        self.setVectorLength(vectorLength)

    def __str__(self):
        return "Angle class with an angle from horizontal of {0} and azimuth of {1}".format(self.angleFromHorizontal,
                                                                                            self.azimuthAngle)

    def setAngleFromHorizontal(self, angleFromHorizontal):
        self.angleFromHorizontal = angleFromHorizontal
        self.updateComponents()

    def setAzimuth(self, azimuthAngle):
        self.azimuthAngle = azimuthAngle
        self.updateComponents()

    def setAngles(self, angleFromHorizontal, azimuthAngle):
        self.angleFromHorizontal = angleFromHorizontal
        self.azimuthAngle = azimuthAngle
        self.updateComponents()

    def addAngle(self, angleObject):
        newAngleXPart = self.xComponent * self.vectorLength + angleObject.xComponent * angleObject.vectorLength
        newAngleYPart = self.yComponent * self.vectorLength + angleObject.yComponent * angleObject.vectorLength
        newAngleZPart = self.zComponent * self.vectorLength + angleObject.zComponent * angleObject.vectorLength

        self.setVectorLength((newAngleXPart ** 2 + newAngleYPart ** 2 + newAngleZPart ** 2) ** 0.5)
        if newAngleXPart == 0 and newAngleYPart == 0:
            if newAngleZPart < 0:
                self.setAngles(-90, 0)
            else:
                self.setAngles(90, 0)
        else:
            # means that length is not 0, so normalize components
            newAngleXPart /= self.vectorLength
            newAngleYPart /= self.vectorLength
            newAngleZPart /= self.vectorLength
            elevationAngle = radiansToDegrees(math.asin(newAngleZPart))
            angleFromEast = radiansToDegrees(
                math.acos(newAngleXPart / (newAngleXPart ** 2 + newAngleYPart ** 2) ** 0.5))
            if newAngleYPart < 0:
                # means that is on bottom part of unit circle, so add 180 degrees
                angleFromEast += 180
            azimuthAngle = 90 - angleFromEast
            if azimuthAngle >= 360:
                azimuthAngle -= 360
            elif azimuthAngle < 0:
                azimuthAngle += 360
            self.setAngles(elevationAngle, azimuthAngle)

    def setVectorLength(self, vectorLength):
        self.vectorLength = vectorLength

    def updateComponents(self):
        elevationAngle = degreesToRadians(self.angleFromHorizontal)  # 0 at horizontal, pi/2 at straight vertical
        angleFromEast = degreesToRadians(
            90 - self.azimuthAngle)  # change from Azimuth (where North is 0, positive going East), to angle (East is 0, positive going North)

        self.xComponent = math.cos(elevationAngle) * math.cos(angleFromEast)
        self.yComponent = math.cos(elevationAngle) * math.sin(angleFromEast)
        self.zComponent = math.sin(elevationAngle)

    def compareToOtherAngle(self, angleObject):
        # how the two angles compare - essentially the dot product (in range of -1 to 1)
        return self.xComponent * angleObject.xComponent + self.yComponent * angleObject.yComponent + self.zComponent * angleObject.zComponent

    def perpendicularToOtherAngle(self, angleObject):
        # how perpendicular this angle is to another angle. Essentially the magnitude of the cross product (in range of 0 to 1)
        # find the cross product
        crossProductX = self.yComponent * angleObject.zComponent - self.zComponent * angleObject.yComponent
        crossProductY = -(self.xComponent * angleObject.zComponent - self.zComponent * angleObject.xComponent)
        crossProductZ = self.xComponent * angleObject.yComponent - self.yComponent * angleObject.xComponent
        crossProductMagnitude = (crossProductX ** 2 + crossProductY ** 2 + crossProductZ ** 2) ** 0.5
        return crossProductMagnitude

class SunlightHandler:
    # the sunlight intensity is in W/m^2, and is the maximum achievable sunlight intensity on a flat surface in the area
    # the sunAngleHandler is an object which contains the following information:
    #   the degree angle from horizontal is the angle of the sun from the horizon (in degrees), from 0 to 90
    #   the degree azimuth is the angle from North, heading east, from the building to the sun
    #       sun is North = 0 degrees, East = 90 degrees, South = 180 degrees, West = 270 degrees
    def __init__(self, maximumSunlightIntensity=930, latitude=33.75, longitude=-84.388,
                 hoursFromGMT=-5, cloudPercentage=None):  # default values are for Atlanta
        self.directionToSun = AngleHandler()
        self.latitude = latitude
        self.longitude = longitude  # longitude is the longitude of the location (range is -180 to 180) (east is positive, west is negative)
        self.hoursFromGMT = hoursFromGMT
        self.normalDirection = AngleHandler(90, 0)
        self.setMaximumSunlightIntensity(maximumSunlightIntensity)
        self.currentDateTime = datetime(2018, 1, 1)
        if cloudPercentage==None:
            self.cloudPercentage = CloudCover()
        else:
            self.cloudPercentage = cloudPercentage

    def __str__(self):
        datetimeStr = "{}/{}/{} {}:{:0>2d}".format(self.currentDateTime.month, self.currentDateTime.day,
                                                   self.currentDateTime.year, self.currentDateTime.hour,
                                                   self.currentDateTime.minute)
        return "({0}) Sun has an intensity of {1} BTU/(hr*ft^2), is {2} degrees above the horizontal, and is {3} degrees east of North.".format(
            datetimeStr, round(self.getSunlightIntensity(), 1), round(self.directionToSun.angleFromHorizontal, 2),
            round(self.directionToSun.azimuthAngle, 2))

    def setSunAngle(self, degreesAboveHorizontal, azimuthAngle):
        # degreesAboveHorizontal is the angle (in degrees) above horizontal of the sun FROM THE VIEWER'S POINT OF VIEW
        # azimuthAngle is the azimuth angle (in degrees) of the sun FROM THE VIEWER'S POINT OF VIEW - North = 0, East = 90, etc.
        self.directionToSun.setAngles(degreesAboveHorizontal, azimuthAngle)

    def getSunlightIntensity(self, faceAngleHandler=None):
        # faceAngleHandler is an AngleHandler instance, which represents the direction of the normal of the wall face
        if faceAngleHandler == None:
            faceAngleHandler = self.normalDirection
        if self.directionToSun.angleFromHorizontal < 0:
            return 0  # means that the sun isn't above the horizon yet
        angleCoefficient = self.directionToSun.compareToOtherAngle(
            faceAngleHandler)  # find the coefficient of sunlight hitting the surface
        if angleCoefficient < 0:
            wallSunHeat = 0  # means that is coming from behind the face, so no sunlight hits the face
        else:
            sunlightIntensityFactor = self.directionToSun.compareToOtherAngle(
                self.normalDirection) * self.cloudPercentage.convertToSunlightIntensityFactor()  # find the coefficient of the intensity of the sun - really intense when it is directly overhead
            wallSunHeat = self.maximumSunlightIntensity * angleCoefficient * sunlightIntensityFactor
        return wallSunHeat

    def setMaximumSunlightIntensity(self, maximumSunlightIntensity):
        # The given value is in W/m^2, and is the max for the specified region.
        # Calculate the angle which gives the maximum intenstiy, and then backtrack into what true solar intensity value will give this maximum at the prescribed angle
        if abs(self.latitude) <= 22.5:
            trueSolarIntensity = maximumSunlightIntensity
        else:
            earthNormal = AngleHandler(abs(self.latitude))
            maxSunPosition = AngleHandler(22.5)
            trueSolarIntensity = maximumSunlightIntensity / earthNormal.compareToOtherAngle(maxSunPosition)

        # After backed into what the true solar value is, convert to BTU / (hr * ft^2)
        self.maximumSunlightIntensity = trueSolarIntensity * 3.412 / 10.7639

    def getSunAngles(self, localTime):
        # given a time in a standard time, this function calculates the angle that the sun should be in the sky
        # calculations supported by http://www.susdesign.com/sunposition/index.php, in the sunposition.js file
        degreesToRadiansConversion = math.pi / 180.
        radiansToDegreeConversion = 180. / math.pi
        if (localTime.month > 2):
            correctedYear = 2009
            correctedMonth = localTime.month - 3
        else:
            correctedYear = 2008
            correctedMonth = localTime.month + 9

        currentMinute = localTime.hour * 60 + localTime.minute
        UT = (currentMinute / 60) - self.hoursFromGMT
        t = ((UT / 24.0) + localTime.day + int(30.6 * correctedMonth + 0.5) + int(
            365.25 * (correctedYear - 1976)) - 8707.5) / 36525.0
        G = self.normalizeAngleTo360(357.528 + 35999.05 * t)
        C = (1.915 * math.sin(G * degreesToRadiansConversion)) + (
                    0.020 * math.sin(2.0 * G * degreesToRadiansConversion))
        L = self.normalizeAngleTo360(280.460 + (36000.770 * t) + C)
        alpha = L - 2.466 * math.sin(2.0 * L * degreesToRadiansConversion) + 0.053 * math.sin(
            4.0 * L * degreesToRadiansConversion)
        obliquity = 23.4393 - 0.013 * t
        declination = math.atan(math.tan(obliquity * degreesToRadiansConversion) * math.sin(
            alpha * degreesToRadiansConversion)) * radiansToDegreeConversion
        eotAdjustment = (L - C - alpha) / 15.0 * 60.0
        solarTimeMinutes = currentMinute - 60 * (
                    self.hoursFromGMT - self.longitude / 15.0) + eotAdjustment  # in minutes; doesn't include DST
        hourAngle = (solarTimeMinutes - 12 * 60) / 4
        altitudeAngle = radiansToDegreeConversion * math.asin((math.cos(
            self.latitude * degreesToRadiansConversion) * math.cos(declination * degreesToRadiansConversion) * math.cos(
            hourAngle * degreesToRadiansConversion)) + (math.sin(self.latitude * degreesToRadiansConversion) * math.sin(
            declination * degreesToRadiansConversion)))
        azimuthMeat = ((math.sin(altitudeAngle * degreesToRadiansConversion) * math.sin(
            self.latitude * degreesToRadiansConversion)) - math.sin(declination * degreesToRadiansConversion)) / (
                                  math.cos(altitudeAngle * degreesToRadiansConversion) * math.cos(
                              self.latitude * degreesToRadiansConversion))
        if azimuthMeat > 1:
            azimuthMeat = 1
        elif azimuthMeat < -1:
            azimuthMeat = -1
        azimuthAngle = 360 - self.normalizeAngleTo360((radiansToDegreeConversion * math.acos(azimuthMeat)) + 180)
        self.directionToSun.setAngles(altitudeAngle, azimuthAngle)
        self.currentDateTime = localTime

    def normalizeAngleTo360(self, angle):
        # constrain the angle between -360 and 360
        if abs(angle) >= 360:
            constrainedAngle = math.fmod(angle, 360)
        else:
            constrainedAngle = angle

        # convert the constrained angle to 0 to 360, if required
        if constrainedAngle < 0:
            return constrainedAngle + 360
        else:
            return constrainedAngle

class CloudCover:
    def __init__(self):
        self.nextStatusCheckHour = 0
        self.nextCloudPercentageCheckHour = 0
        self.currentHour = 0
        self.monthlyChart = [[.42, .1, .1, .07, .31], [.4, .1, .1, .09, .31], [.38, .1, .1, .09, .33],
                        [.33, .1, .1, .11, .36], [.3, .12, .1, .13, .35], [.27, .16, .12, .16, .29],
                        [.24, .20, .15, .14, .27], [.21, .18, .09, .19, .33], [.23, .11, .08, .16, .42],
                        [.26, .08, .07, .11, .48], [.33, .09, .09, .08, .41], [.39, .12, .10, .07, .32]]
        self.percentageOfSkyCovered = [[.7, 1], [.5, .7], [.25, .5], [.05, .25], [0, .05]] # cloudy, mostly cloudy, partly cloudy, mostly clear, clear
        self.defaultCloudCoverRange = (1, 1) # modulated within the increment time
        self.cloudCoverPercentage = 0.

    def resetToTime(self, hourOfDay=0., monthOfSimulation=6, resetCurrentBounds=True):
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        self.currentHour = hourOfDay
        if resetCurrentBounds:
            self.nextStatusCheckHour = hourOfDay - 0.01 # subtract off a bit, so when recalculation occurs, will set the correct value
            self.nextCloudPercentageCheckHour = hourOfDay - 0.01
        self.month = monthOfSimulation

    def incrementTime(self, hourDifference=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        # hourDifference is the amount that the time is changing, in hours
        # update the current time, rolling over if necessary
        self.currentHour += hourDifference
        if self.currentHour >= 24:
            self.currentHour = self.currentHour % 24
            # update the last status checks and last schedule checks here, to avoid rolling over later, and checking at invalid times
            if self.nextStatusCheckHour >= 24:
                self.nextStatusCheckHour = self.nextStatusCheckHour % 24
            if self.nextCloudPercentageCheckHour >= 23.98:
                self.nextCloudPercentageCheckHour = (self.nextCloudPercentageCheckHour % 23) - 0.99

        if self.currentHour > self.nextStatusCheckHour:  # reassign the status
            self.nextStatusCheckHour += round(
                selectRandomFromRangeList(((1, 3), (3, 6), (6, 9), (9, 12)), (0.1, 0.5, 0.3, 0.1)), 3)
            self.defaultCloudCoverRange = selectRandomElement(self.percentageOfSkyCovered, self.monthlyChart[self.month - 1])

        if self.nextCloudPercentageCheckHour < self.currentHour:
            self.percentageOfSkyCovered.remove(self.defaultCloudCoverRange)
            self.outOfRangeCloudRange = selectRandomElement(self.percentageOfSkyCovered)    # choose a different range to be the cloud cover range 2% of the time
            self.cloudCoverPercentage = selectRandomFromRangeList(((self.defaultCloudCoverRange), self.outOfRangeCloudRange), (0.98, 0.02))
            self.percentageOfSkyCovered.append(self.defaultCloudCoverRange)

        if self.cloudCoverPercentage > self.defaultCloudCoverRange[1] or self.cloudCoverPercentage < self.defaultCloudCoverRange[0]:
            if self.currentHour > self.nextCloudPercentageCheckHour:
                self.nextCloudPercentageCheckHour += round(selectRandomFromRangeList(((0., 0.05), (0.05, 0.2), (0.2, 0.3), (0.3, 0.5)), (0.3, 0.45, 0.15, 0.1)), 3)
                self.cloudCoverPercentage = selectRandomFromRange(self.outOfRangeCloudRange[0], self.outOfRangeCloudRange[1])
            else:
                self.cloudCoverPercentage = selectRandomFromRange(self.outOfRangeCloudRange[0], self.outOfRangeCloudRange[1])
        else:
            self.nextCloudPercentageCheckHour = self.currentHour - 0.01
        self.convertToSunlightIntensityFactor()

    def convertToSunlightIntensityFactor(self):
        return 1. - 0.5 * self.cloudCoverPercentage

class OutsideTemperature:
    def __init__(self):
        self.nextStatusCheckHour = 0
        self.currentHour = 0
        self.monthlyLowTempChart = [[[61., 54.], [54., 46.], [46., 42.9], [42.9, 34.6], [34.6, 25.9], [25.9, 29.], [29., 17.], [17., 6.]],
                                 [[65., 53.], [53., 44.], [44., 48.9], [48.9, 38.1], [38.1, 30.8], [30.8, 32.], [32., 24.], [24., 15.]],
                                 [[65., 60.], [60., 52.], [52., 54.7], [54.7, 45.4], [45.4, 39.4], [39.4, 35.], [35., 29.], [29., 21.]],
                                 [[69., 64.], [64., 58.], [58., 58.], [58., 52.8], [52.8, 47.4], [47.4, 45.], [45., 37.], [37., 28.]],
                                 [[73., 70.], [70., 66.], [66., 65.9], [65.9, 61.4], [61.4, 56.8], [56.8, 54.], [54., 48.], [48., 41.]],
                                 [[77., 74.], [74., 70.], [70., 72.3], [72.3, 68.9], [68.9, 66], [66., 67.], [67., 61.], [61., 53.]],
                                 [[79., 76.], [76., 73.], [73., 73.5], [73.5, 71.4], [71.4, 69.3], [69.3, 70.], [70., 65.], [65., 59.]],
                                 [[82., 76.], [76., 73.], [73., 75.2], [75.2, 71.4], [71.4, 68.2], [68.2, 72.], [72., 65.], [65., 56.]],
                                 [[75., 73.], [73., 71.], [71., 69.3], [69.3, 65.6], [65.6, 61.8], [61.8, 62.], [62., 53.], [53., 43.]],
                                 [[73., 67.], [67., 61.], [61., 59.4], [59.4, 54.5], [54.5, 49.5], [49.5, 45.], [45., 39.], [39., 32.]],
                                 [[69., 60.], [60., 51.], [51., 49.5], [49.5, 44.3], [44.3, 38.1], [38.1, 36.], [36., 29.], [29., 21.]],
                                 [[66., 54.], [54., 42.], [42., 49.5], [49.5, 38.], [38., 29.], [29., 30.], [30., 24.], [24., 13.]]]
        self.monthlyHighTempChart = [[[78., 71.], [71., 66.], [66., 61.2], [61.2, 53.3], [53.3, 47.3], [47.3, 45.], [45., 33.], [33., 26.]],
                                     [[79., 72.], [72., 64.], [64., 66.5], [66.5, 57.6], [57.6, 48.5], [48.5, 50.], [50., 39.], [39., 28.]],
                                     [[87., 80.], [80., 72.], [72., 74.2], [74.2, 65.5], [65.5, 58.9], [58.9, 53.], [53., 46.], [46., 35.]],
                                     [[88., 85.], [85., 78.], [78., 77.7], [77.7, 73.3], [73.3, 69.8], [69.8, 64.], [64., 56.], [56., 46.]],
                                     [[93., 90.], [90., 86.], [86., 83.8], [83.8, 80.6], [80.6, 76.7], [76.7, 77.], [77., 67.], [67., 56.]],
                                     [[106., 94.], [94., 88.], [88., 91.9], [91.9, 87.3], [87.3, 82.6], [82.6, 86.], [86., 78.], [78., 64.]],
                                     [[105., 95.], [95., 90.], [90., 94.1], [94.1, 89.1], [89.1, 84.7], [84.7, 89.], [89., 80.], [80., 73.]],
                                     [[104., 96.], [96., 92.], [92., 96], [96., 88.7], [88.7, 85.], [85., 86.], [86., 79.], [79., 66.]],
                                     [[96., 91.], [91., 85.], [85., 88.9], [88.9, 82.9], [82.9, 78.7], [78.7, 77.], [77., 70.], [70., 64.]],
                                     [[88., 85.], [85., 81.], [81., 81.1], [81.1, 74.], [74., 69.3], [69.3, 68.], [68., 57], [57., 45.]],
                                     [[84., 78.], [78., 74.], [74., 70.8], [70.8, 64.3], [64.3, 59.1], [59.1, 59.], [59., 47.], [47., 37.]],
                                     [[77., 71.], [71., 64.], [64., 65.7], [65.7, 55.8], [55.8, 45.4], [45.4, 50.], [50., 40.], [40., 29.]]]
        self.monthlyPercentLikelyhood = [0.05, 0.1, 0.15, 0.2, 0.2, 0.15, 0.1, 0.05]
        self.currentTemp = 0

    def resetToTime(self, hourOfDay=0., monthOfSimulation=6, resetCurrentBounds=True):
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        self.currentHour = hourOfDay
        self.month = monthOfSimulation
        if resetCurrentBounds:
            self.lowTemp = selectRandomFromRangeList(self.monthlyLowTempChart[self.month - 1], self.monthlyPercentLikelyhood)
            self.lowTempTime = selectRandomFromRangeList(((0, 2), (2, 4), (4, 6), (6, 8)), (0.05, 0.1, 0.35, 0.5))
            self.highTemp = selectRandomFromRangeList(self.monthlyHighTempChart[self.month - 1], self.monthlyPercentLikelyhood)
            self.highTempTime = selectRandomFromRangeList(((12, 14), (14, 16), (16, 18), (18, 20)), (0.05, 0.35, 0.5, 0.1))
            self.timeDiff = self.highTempTime - self.lowTempTime
            self.oldLowTemp = self.lowTemp
            self.oldHighTemp = self.highTemp
            if self.currentHour > self.lowTempTime and self.currentHour < self.highTempTime:
                self.nextStatusCheckHour = self.highTempTime
            else:
                self.nextStatusCheckHour = self.lowTempTime

    def incrementTime(self, hourDifference=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        # hourDifference is the amount that the time is changing, in hours
        # update the current time, rolling over if necessary
        self.currentHour += hourDifference
        if self.currentHour >= 24:
            self.currentHour = self.currentHour % 24
            # update the last status checks and last schedule checks here, to avoid rolling over later, and checking at invalid times
            if self.nextStatusCheckHour >= 24:
                self.nextStatusCheckHour = self.nextStatusCheckHour % 24

        if self.currentHour > self.nextStatusCheckHour:  # reassign the status
            if self.currentHour >= 0. and self.currentHour <= 8.:
                self.oldLowTemp = self.lowTemp
                self.timeDiff = self.highTempTime - self.lowTempTime
                self.lowTemp = selectRandomFromRangeList(self.monthlyLowTempChart[self.month - 1],
                                                         self.monthlyPercentLikelyhood)
                self.lowTempTime = selectRandomFromRangeList(((0, 2), (2, 4), (4, 6), (6, 8)), (0.05, 0.1, 0.35, 0.5))
                self.nextStatusCheckHour = self.highTempTime             # reset next checkpoint and low temperature
            elif self.currentHour >= 12. and self.nextStatusCheckHour != self.lowTempTime:
                self.oldHighTemp = self.highTemp
                self.timeDiff = (24. - self.highTempTime) + self.lowTempTime
                self.highTemp = selectRandomFromRangeList(self.monthlyHighTempChart[self.month - 1],
                                                          self.monthlyPercentLikelyhood)
                self.highTempTime = selectRandomFromRangeList(((12, 14), (14, 16), (16, 18), (18, 20)),
                                                              (0.05, 0.35, 0.5, 0.1))
                self.nextStatusCheckHour = self.lowTempTime              # reset next checkpoint and high temperature

        if self.nextStatusCheckHour == self.highTempTime:
            self.currentTemp = (self.highTemp - self.oldLowTemp) * math.sin((1. - ((self.nextStatusCheckHour - self.currentHour) / self.timeDiff)) * (math.pi / 2)) + self.oldLowTemp
        elif self.nextStatusCheckHour == self.lowTempTime:
            if self.currentHour >= 12. and self.currentHour <= 24.:
                self.currentTemp = (self.lowTemp - self.oldHighTemp) * math.sin(
                    (1. - ((self.nextStatusCheckHour + (24. - self.currentHour)) / self.timeDiff)) * (math.pi / 2)) + self.oldHighTemp
            elif self.currentHour >= 0. and self.currentHour <= 12.:
                self.currentTemp = (self.lowTemp - self.oldHighTemp) * math.sin(
                    (1. - ((self.nextStatusCheckHour - self.currentHour) / self.timeDiff)) * (
                                math.pi / 2)) + self.oldHighTemp
        return self.currentTemp

class OutsideHumidity:
    def __init__(self):
        self.nextStatusCheckHour = 0
        self.currentHour = 0
        self.monthlyLowHumidityChart = [[[68.0, 65.0], [65.0, 62.0], [62.0, 59.0], [59.0, 56.0], [56.0, 53.0], [53.0, 50.0]],
                                         [[64.0, 61.0], [61.0, 58.0], [58.0, 55.0], [55.0, 52.0], [52.0, 49.0], [49.0, 46.0]],
                                         [[61.0, 58.0], [58.0, 55.0], [55.0, 52.0], [52.0, 49.0], [49.0, 46.0], [46.0, 43.0]],
                                         [[59.0, 56.0], [56.0, 53.0], [53.0, 50.0], [50.0, 47.0], [47.0, 44.0], [44.0, 41.0]],
                                         [[63.0, 60.0], [60.0, 57.0], [57.0, 54.0], [54.0, 51.0], [51.0, 48.0], [48.0, 45.0]],
                                         [[66.0, 63.0], [63.0, 60.0], [60.0, 57.0], [57.0, 54.0], [54.0, 51.0], [51.0, 48.0]],
                                         [[68.0, 65.0], [65.0, 62.0], [62.0, 59.0], [59.0, 56.0], [56.0, 53.0], [53.0, 50.0]],
                                         [[69.0, 66.0], [66.0, 63.0], [63.0, 60.0], [60.0, 57.0], [57.0, 54.0], [54.0, 51.0]],
                                         [[68.0, 65.0], [65.0, 62.0], [62.0, 59.0], [59.0, 56.0], [56.0, 53.0], [53.0, 50.0]],
                                         [[63.0, 60.0], [60.0, 57.0], [57.0, 54.0], [54.0, 51.0], [51.0, 48.0], [48.0, 45.0]],
                                         [[65.0, 62.0], [62.0, 59.0], [59.0, 56.0], [56.0, 53.0], [53.0, 50.0], [50.0, 47.0]],
                                         [[68.0, 65.0], [65.0, 62.0], [62.0, 59.0], [59.0, 56.0], [56.0, 53.0], [53.0, 50.0]]]

        self.monthlyHighHumidityChart = [[[87.0, 84.0], [84.0, 81.0], [81.0, 78.0], [78.0, 75.0], [75.0, 72.0], [72.0, 69.0]],
                                          [[86.0, 83.0], [83.0, 80.0], [80.0, 77.0], [77.0, 74.0], [74.0, 71.0], [71.0, 68.0]],
                                          [[87.0, 84.0], [84.0, 81.0], [81.0, 78.0], [78.0, 75.0], [75.0, 72.0], [72.0, 69.0]],
                                          [[88.0, 85.0], [85.0, 82.0], [82.0, 79.0], [79.0, 76.0], [76.0, 73.0], [73.0, 70.0]],
                                          [[91.0, 88.0], [88.0, 85.0], [85.0, 82.0], [82.0, 79.0], [79.0, 76.0], [76.0, 73.0]],
                                          [[93.0, 90.0], [90.0, 87.0], [87.0, 84.0], [84.0, 81.0], [81.0, 78.0], [78.0, 75.0]],
                                          [[97.0, 94.0], [94.0, 91.0], [91.0, 88.0], [88.0, 85.0], [85.0, 82.0], [82.0, 79.0]],
                                          [[98.0, 95.0], [95.0, 92.0], [92.0, 89.0], [89.0, 86.0], [86.0, 83.0], [83.0, 80.0]],
                                          [[96.0, 93.0], [93.0, 90.0], [90.0, 87.0], [87.0, 84.0], [84.0, 81.0], [81.0, 78.0]],
                                          [[93.0, 90.0], [90.0, 87.0], [87.0, 84.0], [84.0, 81.0], [81.0, 78.0], [78.0, 75.0]],
                                          [[88.0, 85.0], [85.0, 82.0], [82.0, 79.0], [79.0, 76.0], [76.0, 73.0], [73.0, 70.0]],
                                          [[86.0, 83.0], [83.0, 80.0], [80.0, 77.0], [77.0, 74.0], [74.0, 71.0], [71.0, 68.0]]]

        self.monthlyPercentLikelyhood = [0.05, 0.1, 0.35, 0.35, 0.1, 0.05]

    def resetToTime(self, hourOfDay=0., monthOfSimulation=6, resetCurrentBounds=True):
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        self.currentHour = hourOfDay
        self.month = monthOfSimulation
        if resetCurrentBounds:
            self.lowHumidity = selectRandomFromRangeList(self.monthlyLowHumidityChart[self.month - 1], self.monthlyPercentLikelyhood)
            self.lowHumidityTime = selectRandomFromRangeList(((13, 15), (15, 17), (17, 19)), (0.2, 0.6, 0.2))
            self.highHumidity = selectRandomFromRangeList(self.monthlyHighHumidityChart[self.month - 1], self.monthlyPercentLikelyhood)
            self.highHumidityTime = selectRandomFromRangeList(((4, 6), (6, 8), (8, 10)), (0.2, 0.6, 0.2))
            self.timeDiff = self.highHumidityTime - self.lowHumidityTime
            self.oldLowHumidity = self.lowHumidity
            self.oldHighHumidity = self.highHumidity
            if self.currentHour > self.highHumidityTime and self.currentHour < self.lowHumidityTime:
                self.nextStatusCheckHour = self.lowHumidityTime
            else:
                self.nextStatusCheckHour = self.highHumidityTime

    def incrementTime(self, hourDifference=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        # hourDifference is the amount that the time is changing, in hours
        # update the current time, rolling over if necessary
        self.currentHour += hourDifference
        if self.currentHour >= 24:
            self.currentHour = self.currentHour % 24
            # update the last status checks and last schedule checks here, to avoid rolling over later, and checking at invalid times
            if self.nextStatusCheckHour >= 24:
                self.nextStatusCheckHour = self.nextStatusCheckHour % 24

        if self.currentHour > self.nextStatusCheckHour:  # reassign the status
            if self.currentHour >= 4. and self.currentHour <= 10.:
                self.oldHighHumidity = self.highHumidity
                self.timeDiff = self.lowHumidityTime - self.highHumidityTime
                self.highHumidity = selectRandomFromRangeList(self.monthlyHighHumidityChart[self.month - 1],
                                                         self.monthlyPercentLikelyhood)
                self.highHumidityTime = selectRandomFromRangeList(((4, 6), (6, 8), (8, 10)), (0.2, 0.6, 0.2))
                self.nextStatusCheckHour = self.lowHumidityTime             # reset next checkpoint and low temperature
            elif self.currentHour >= 13. and self.nextStatusCheckHour != self.highHumidityTime:
                self.oldLowHumidity = self.lowHumidity
                self.timeDiff = (24. - self.lowHumidityTime) + self.highHumidityTime
                self.lowHumidity = selectRandomFromRangeList(self.monthlyLowHumidityChart[self.month - 1],
                                                          self.monthlyPercentLikelyhood)
                self.lowHumidityTime = selectRandomFromRangeList(((13, 15), (15, 17), (17, 19)), (0.2, 0.6, 0.2))
                self.nextStatusCheckHour = self.highHumidityTime             # reset next checkpoint and high temperature

        if self.nextStatusCheckHour == self.lowHumidityTime:
            self.currentHumidity = (self.lowHumidity - self.oldHighHumidity) * math.sin((1. - ((self.nextStatusCheckHour - self.currentHour) / self.timeDiff)) * (math.pi / 2)) + self.oldHighHumidity
        elif self.nextStatusCheckHour == self.highHumidityTime:
            if self.currentHour >= 13. and self.currentHour <= 24.:
                self.currentHumidity = (self.highHumidity - self.oldLowHumidity) * math.sin(
                    (1. - ((self.nextStatusCheckHour + (24. - self.currentHour)) / self.timeDiff)) * (math.pi / 2)) + self.oldLowHumidity
            elif self.currentHour >= 0. and self.currentHour <= 10.:
                self.currentHumidity = (self.highHumidity - self.oldLowHumidity) * math.sin(
                    (1. - ((self.nextStatusCheckHour - self.currentHour) / self.timeDiff)) * (
                                math.pi / 2)) + self.oldLowHumidity
        return self.currentHumidity

class WindHandler:
    # class which describes the wind's effect on a building
    # wind speed as parameter is in miles per hour, but is stored internally as meters per second
    def __init__(self, windSpeed=0, windAzimuth=0):
        self.windDirection = AngleHandler(0)  # set the angle above horizontal to 0
        self.setWindCharacteristics(windSpeed, windAzimuth)

    def setWindCharacteristics(self, windSpeed=None, windAzimuth=None):
        # wind speed is in miles per hour, wind azimuth is in degrees from North, going East (North = 0, East = 90, etc)
        if windSpeed != None:
            self.windSpeed = windSpeed * 0.44704  # convert from miles per hour to meters per second
            self.windDirection.setVectorLength(self.windSpeed)
        if windAzimuth != None:
            self.windDirection.setAzimuth(windAzimuth)

    def getWindConvectionCoefficient(self, faceAngleHandler):
        # faceAngleHandler is an AngleHandler instance, which represents the direction of the normal of the wall face
        # convection coefficient will be when the wind is blowing across the face, which means wind direction is perpindicular to face normal
        # First, find the dissimilarity (perpendicularity) of the two angles
        relativePerpendicularWindspeed = self.windSpeed * self.windDirection.perpendicularToOtherAngle(faceAngleHandler)

        # second, calculate the convection coefficient factor, in W/(m^2 * K)
        if relativePerpendicularWindspeed > 2:
            # if windspeed is between 2 m/s and 20 m/s, the following equation can be used:
            # equation from https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
            convectiveCoefficient = 10.45 - relativePerpendicularWindspeed + 10 * relativePerpendicularWindspeed ** 0.5
        else:
            # interpolate between 5 and 22.5 (the lowerbound limit of the above equation
            lowerCoefficient = 5
            lowerWindspeed = 0
            upperCoefficient = 22.5
            upperWindspeed = 2
            convectiveCoefficient = lowerCoefficient + (upperCoefficient - lowerCoefficient) * (
                        relativePerpendicularWindspeed - lowerWindspeed) / (upperWindspeed - lowerWindspeed)
        return convectiveCoefficient * 0.1761  # convert from W / (m^2 * degC) to BTU / (hr * ft^2 * degF)

    def addWindDirection(self, additionalWindSpeed=0, additionalWindAzimuth=0):
        additionalWindAngle = AngleHandler(0, additionalWindAzimuth, additionalWindSpeed)
        self.windDirection.addAngle(additionalWindAngle)

    def randomlyIncrementWind(self, incrementMinutes=5.):
        maxDifference = max(incrementMinutes / 5 * 2, 5)
        additionalWindSpeed = selectRandomFromRange(0, maxDifference)
        additionalWindAzimuth = selectRandomFromRange(0, 360)
        self.addWindDirection(additionalWindSpeed, additionalWindAzimuth)

    def randomlySetWind(self):
        windSpeed = round(
            selectRandomFromRangeList(((0, 2), (2, 5), (5, 10), (10, 20), (20, 40)), (0.2, 0.2, 0.3, 0.25, 0.05)), 2)
        windAzimuth = selectRandomFromRange(0, 360)
        self.setWindCharacteristics(windSpeed, windAzimuth)

class OutsideConditionsHandler:
    def __init__(self):
        self.outdoorTemperature = 0  # in degF
        self.outdoorHumidity = 70  # in % relative humidity (0 to 100)
        self.CloudCover = CloudCover()
        self.OutsideTemp = OutsideTemperature()
        self.OutsideHumidity = OutsideHumidity()
        self.sunlightData = SunlightHandler(cloudPercentage=self.CloudCover)
        self.windData = WindHandler()
        self.windData.randomlySetWind()

    def setDesignHeatingConditions(self):
        self.outdoorTemperature = 20
        self.outdoorHumidity = 60
        self.windData.setWindCharacteristics(windSpeed=10)
        self.setTime(2018, 1, 15, 0, 0)  # Midnight in the middle of January

    def setDesignCoolingConditions(self):
        self.outdoorTemperature = 99
        self.outdoorHumidity = 80
        self.windData.setWindCharacteristics(windSpeed=10)
        self.setTime(2018, 7, 15, 14, 0)  # Afternoon (2P) in middle of July

    def setTime(self, currYear=2018, currMonth=4, currDay=1, currHour=12, currMinute=0, changeBounds = True):
        self.modelTime = datetime(currYear, currMonth, currDay, currHour, currMinute)
        self.sunlightData.getSunAngles(self.modelTime)
        self.windData.randomlySetWind()
        self.CloudCover.resetToTime(currHour, currMonth, changeBounds)
        self.OutsideTemp.resetToTime(currHour, currMonth, changeBounds)
        self.OutsideHumidity.resetToTime(currHour, currMonth, changeBounds)

    def incrementTime(self, numMinutes=1):
        #Simulation runs every minute, data log operates every 15 minutes (completely separate)
        self.modelTime += timedelta(minutes=numMinutes)
        self.sunlightData.getSunAngles(self.modelTime)
        self.windData.randomlyIncrementWind(numMinutes)
        self.CloudCover.incrementTime(numMinutes)
        self.outdoorTemperature = self.OutsideTemp.incrementTime(numMinutes)
        self.outdoorHumidity = self.OutsideHumidity.incrementTime(numMinutes)

class RoomFace:
    # class which represents a wall in a room
    def __init__(self, faceArea,  # the area of the room face
                 constructionMaterialParameters,
                 # instance of ConstructionMaterialValues which contains construction information
                 degreeAngleFromHorizontal=0,  # normal direction of wall in degrees from horizontal (0 - 90 degrees)
                 degreeAngleAzimuth=0,
                 # normal direction of wall, in azimuth (north is 0, east is 90, etc), pointing AWAY from room
                 isExterior=False,
                 faceName="Wall Face"):  # name of the face
        self.constructionMaterialParameters = constructionMaterialParameters
        self.faceArea = faceArea
        self.faceName = faceName
        self.isExterior = (isExterior == True)
        self.faceNormalHandler = AngleHandler(degreeAngleFromHorizontal,
                                              degreeAngleAzimuth)  # azimuth angle from north (going East) of wall's normal direction, in degrees
        self.generateWallParameters()
        self.exteriorTemp = 70

    def generateWallParameters(self):
        if self.isExterior:
            if self.faceNormalHandler.angleFromHorizontal > 45:  # means that is probably roof
                self.percentageWindow = 0.
                wallInsulation = self.constructionMaterialParameters.roofRValue
            elif self.faceNormalHandler.angleFromHorizontal < -45:  # means that is probably a floor
                self.percentageWindow = 0.
                wallInsulation = self.constructionMaterialParameters.floorRValue
            else:
                self.percentageWindow = round(
                    selectRandomFromRangeList(((0, 0), (0.1, 0.25), (0.25, 0.4), (0.4, 0.8), (0.9, 1)),
                                              (0.3, 0.1, 0.3, 0.2, 0.1)), 2)
                wallInsulation = self.constructionMaterialParameters.outerWallRValue
        else:
            self.percentageWindow = 0.
            wallInsulation = self.constructionMaterialParameters.innerWallRValue
        windowInsulation = self.constructionMaterialParameters.windowRValue
        self.windowArea = self.faceArea * self.percentageWindow  # figure out how much of the
        self.wallArea = self.faceArea * (1 - self.percentageWindow)
        self.wallInsulation = wallInsulation
        self.effectiveRValue = 1 / (
                    self.percentageWindow / windowInsulation + (1 - self.percentageWindow) / wallInsulation)

    def getHeatFlowToRoom(self, roomTemperature, otherSideTemperature, outsideConditions):
        # roomTemperature is the temperature within the room
        # otherSideTemperature is the interior temperature of the room on the other side of the wall (floor average)
        # outsideConditions is an instance of WindHandler, which contains information about the outdoor conditions
        if not self.isExterior:
            # interior wall, so determine just the heat flow through the wall
            return (otherSideTemperature - roomTemperature) * self.faceArea / self.effectiveRValue
        else:
            # get the insolation value of the sunlight on this face
            sunlightIntensity = outsideConditions.sunlightData.getSunlightIntensity(
                self.faceNormalHandler)  # in BTU / (hr * ft^2)

            # get the outer wall temperature
            outerWallTemp = self.getTemperatureOfExteriorWall(roomTemperature, sunlightIntensity, outsideConditions)

            # get the heat from sunlight through the window
            heatFromSunThroughWindow = self.windowArea * sunlightIntensity * self.constructionMaterialParameters.windowTransmittance

            # get the conduction through the wall and windows
            wallConductiveHeatFlow = self.wallArea * (outerWallTemp - roomTemperature) / self.wallInsulation
            windowConductiveHeatFlow = self.windowArea * (outsideConditions.outdoorTemperature - roomTemperature) \
                                       / self.constructionMaterialParameters.windowRValue

            # return the sum of all of the heat coming into / going out of the room
            return heatFromSunThroughWindow + wallConductiveHeatFlow + windowConductiveHeatFlow

    def getTemperatureOfExteriorWall(self, roomTemperature, sunlightIntensity, outsideConditions):
        if self.faceNormalHandler.angleFromHorizontal > 45:  # means is facing upwards, so use roof emmittance
            faceEmittance = self.constructionMaterialParameters.roofEmittance
        else:
            faceEmittance = self.constructionMaterialParameters.outerWallEmittance
        convectionCoefficient = outsideConditions.windData.getWindConvectionCoefficient(self.faceNormalHandler)
        outsideTemp = outsideConditions.outdoorTemperature
        self.exteriorTemp = solveForExteriorTempOfWall(sunlightIntensity, faceEmittance, convectionCoefficient,
                                                       self.wallInsulation,
                                                       roomTemperature, outsideTemp, self.exteriorTemp, threshold=0.001)
        return self.exteriorTemp

class OccupancySchedule:
    def __init__(self):
        self.softStartHour = 5
        self.fullOccupancyStartHour = 8
        self.fullOccupancyStopHour = 17
        self.softStopHour = 20
        self.percentAlwaysActive = 0.05  # probability that an area is occupied during the night (5% chance it will be occupied during this time frame)
        self.percentShoulder = 0.3  # probability that an area is occupied during the shoulder periods (5am to 8am and 5pm to 8pm)
        self.percentFullOccupancy = 0.9  # probability that an area is occupied during the full occupancy period

    def getOccupiedStatus(self, hourOfDay):
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        # get the percent likelihood that the area is occupied
        percentLikelihood = self.percentFullOccupancy
        if hourOfDay < self.softStartHour or hourOfDay > self.softStopHour:
            # within the 24 hours period, so return the probability of being 24/7
            percentLikelihood = self.percentAlwaysActive
        elif hourOfDay < self.fullOccupancyStartHour or hourOfDay > self.fullOccupancyStopHour:
            # within the shoulder period, so return the probability of being in occupied during the shoulder
            percentLikelihood = self.percentShoulder

        # randomly select, given the likelihood of occupancy, whether the area is occupied
        return selectRandomElement((True, False), (percentLikelihood, 1. - percentLikelihood))

    def getAirHandlerOccupiedStatus(self, hourofDay):
        if hourofDay >= (self.fullOccupancyStartHour - 1) and hourofDay < self.fullOccupancyStopHour:   # for some reason the hour of the day is non-inclusive, so subract 1 to make it include the start hour
            return 1
        else:
            return 0

class RoomUsageStatus:
    def __init__(self, designLighting, designOccupantHeatLoad, designBaseLoad, occupancySchedule=None):
        if occupancySchedule is None:
            self.occupancySchedule = OccupancySchedule()
        else:
            self.occupancySchedule = occupancySchedule
        self.probabilityOccupiedDuringSchedule = 0.8
        self.nextScheduleCheckHour = 0
        self.scheduleCheckIncrements = 1.  # check every hour
        self.nextStatusCheckHour = 0
        self.currentHour = 0
        self.occupiedFromSchedule = True
        self.isOccupied = True
        self.AHUoccupiedFromSchedule = 1
        self.designAmounts = {"LightingDesign": designLighting, "PeopleDesign": designOccupantHeatLoad,
                              "BaseLoadDesign": designBaseLoad}
        self.status = {"LightingFactor": 0, "PeopleFactor": 0, "BaseLoadFactor": 0}
        self.loads = {"LightingLoad": 0, "PeopleLoad": 0, "BaseLoadLoad": 0}

    def resetToTime(self, hourOfDay=0.):
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        self.currentHour = hourOfDay
        self.nextScheduleCheckHour = hourOfDay - 0.01  # subtract off a bit, so when recalculation occurs, will set the correct value
        self.nextStatusCheckHour = hourOfDay - 0.01

    def incrementTime(self, hourDifference=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        # hourDifference is the amount that the time is changing, in hours
        # update the current time, rolling over if necessary
        self.currentHour += hourDifference
        if self.currentHour >= 24:
            self.currentHour = self.currentHour % 24
            # update the last status checks and last schedule checks here, to avoid rolling over later, and checking at invalid times
            if self.nextScheduleCheckHour >= 24:
                self.nextScheduleCheckHour = self.nextScheduleCheckHour % 24
            if self.nextStatusCheckHour >= 24:
                self.nextStatusCheckHour = self.nextStatusCheckHour % 24

        if self.currentHour > self.nextScheduleCheckHour:  # check the schedule to see if the room is occupied
            self.occupiedFromSchedule = self.occupancySchedule.getOccupiedStatus(self.currentHour)
            self.AHUoccupiedFromSchedule = self.occupancySchedule.getAirHandlerOccupiedStatus(self.currentHour)
            self.nextScheduleCheckHour += self.scheduleCheckIncrements
            self.isOccupied = self.occupiedFromSchedule
            self.getAHUstatus()

        if self.currentHour > self.nextStatusCheckHour:  # reassign the status
            statusDuration = round(
                selectRandomFromRangeList(((5, 15), (15, 30), (30, 60), (60, 180)), (0.1, 0.3, 0.5, 0.1)), 0)
            self.nextStatusCheckHour += statusDuration / 60.
            if self.occupiedFromSchedule:
                # select a random weighting of each of the heat load factors
                self.status["LightingFactor"] = selectRandomFromRangeList(
                    ((0.25, 0.7), (0.7, 0.9), (0.9, 1.1), (1.1, 1.3)), (0.05, 0.8, 0.1, 0.05))
                self.status["PeopleFactor"] = selectRandomFromRangeList(
                    ((0.25, 0.7), (0.7, 0.9), (0.9, 1.1), (1.1, 1.5), (1.5, 3)), (0.25, 0.3, 0.3, 0.1, 0.05))
                self.status["BaseLoadFactor"] = selectRandomFromRangeList(
                    ((0.25, 0.7), (0.7, 0.9), (0.9, 1.1), (1.1, 1.5)), (0.25, 0.35, 0.3, 0.1))
            else:
                self.status["LightingFactor"] = 0.
                self.status["PeopleFactor"] = 0.
                self.status["BaseLoadFactor"] = 0.

            self.loads["LightingLoad"] = round(self.designAmounts["LightingDesign"] * self.status["LightingFactor"],
                                               2)  # lighting load, in kW
            self.loads["PeopleLoad"] = round(self.designAmounts["PeopleDesign"] * self.status["PeopleFactor"],
                                             2)  # occupant heat load, in kW
            self.loads["BaseLoadLoad"] = round(self.designAmounts["BaseLoadDesign"] * self.status["BaseLoadFactor"],
                                               2)  # base load, in kW

    def getAHUstatus(self):
        return self.AHUoccupiedFromSchedule

    def getLightLoad(self):
        return self.loads["LightingLoad"]  # lighting load, in kW

    def getLightingFactor(self):
        return self.status["LightingFactor"]

    def getBaseLoad(self):
        return self.loads["BaseLoadLoad"]

    def calculateLoads(self):
        lightLoad = self.loads["LightingLoad"]  # lighting load, in kW
        occupantLoad = self.loads["PeopleLoad"]  # occupant heat load, in kW
        baseLoad = self.loads["BaseLoadLoad"]  # base load, in kW
        return lightLoad + occupantLoad + baseLoad

class RoomTemperatureSetpoints:
    def __init__(self):
        self.designOccCoolingStpt = 72.0
        self.designOccHeatingStpt = 68.0
        self.designUnoccCoolingStpt = 85.0
        self.designUnoccHeatingStpt = 65.0

        self.occupiedCoolingStpt = roundToNearest(selectRandomFromRange(71, 75), 0.5)  # 76.0
        self.occupiedHeatingStpt = roundToNearest(selectRandomFromRange(66, 70), 0.5)  # 72.0
        self.unoccupiedCoolingStpt = roundToNearest(selectRandomFromRange(90., 80.), 0.5)# 85.0
        self.unoccupiedHeatingStpt = roundToNearest(selectRandomFromRange(70., 60.), 0.5)# 65.0

        self.keepDefaultProbability = 0.5

    def copyRoomSetpoints(self):
        shouldUseDefault = selectRandomElement([True, False], [self.keepDefaultProbability, 1-self.keepDefaultProbability])
        if shouldUseDefault:
            return self
        else:
            return RoomTemperatureSetpoints()


class AirHandler:
    def __init__(self, AHU_Name, building, floor, randomizeSensors=True, setUpWorking=None):
        self.AHU_Name = AHU_Name
        self.design_AHU_Cooling_DAT = selectRandomElement((46, 52, 55), (0.1, 0.1, 0.8))  # AHU design DAT in degF
        self.design_TU_Heating_DAT = 110
        self.design_plenum_temperature = round(selectRandomFromRange(75, 85), 0)  # plenum design temp in degF
        self.design_OA_Percentage = 0.2

        # set up relationships to other objects
        self.building = building
        self.floor = floor
        self.terminalUnitList = []
        self.roomList = []
        self.numberTerminalUnitsNeedCooling = 0
        self.numberTerminalUnitsNeedHeating = 0

        self.outsideConsitions = OutsideConditionsHandler()

        self.issueList = []
        self.pointsList = []
        self.pointNameList = []
        self.setUpCWVOutput(setUpWorking=setUpWorking)
        self.setUpDATSensor(setUpWorking=setUpWorking)
        self.setUpRATSensor(setUpWorking=setUpWorking)
        self.setUpStaticPressureSensor(setUpWorking=setUpWorking)
        self.setUpDischargeAirflowSensor(setUpWorking=setUpWorking)
        self.setUpFanSpeed(setUpWorking=setUpWorking)
        self.setUpOutsideTemperature(setUpWorking=setUpWorking)
        self.setUpOutsideHumidity(setUpWorking=setUpWorking)
        self.coolingOutput = coolingCWVOutput()
        self.fanspeedoutput = AHUFanOutput()
        self.StaticPressureSP = 0.7
        self.setUpStaticPressureSPsensor(setUpWorking=setUpWorking)
        self.DATSP = 55.
        self.setUpAHUDATSP(setUpWorking=setUpWorking)
        self.designDTacrosscoil = 10.
        self.CWstartTemp = 45.
        self.RATdeltaT = 3.
        self.CHWstartDeltaT = 3.
        self.designDuctPressure = 0.7
        self.designDuctPressureBeforePIU = 0.5
        self.designDuctPressureAfterPIU = 0.2
        self.Duct_Static_Pressure_BeforePIU = self.designDuctPressureBeforePIU
        self.Duct_Static_Pressure_AfterPIU = self.designDuctPressureAfterPIU
        self.isOccupied = None
        self.setOccupancy(True)
        self.resetSimulation()
        if randomizeSensors:
            self.randomizeSensorStatus()


    def changeStaticPressureSP(self, newPressure):
        self.StaticPressureSP = newPressure

    def changeDATSP(self, newTemp):
        self.DATSP = newTemp

    def sizeAHU(self):
        designCoolingLoad, designHeatingLoad = self.calculateCoolingHeatingRequirements()
        self.CWV.sizeAHUcoils(designCoolingLoad, designHeatingLoad)
        designAirflow = self.calculateDesignAirflow()
        self.sizeAHUFans(designAirflow)

    def sizeAHUFans(self, designAirflow):
        self.designAirflow = designAirflow
        # calculate the design pressure drop between the static pressure sensor and the AHU
        # assume that the AHU is in the middle of the building, and goes from the middle to a corner (so half of width and half of length)
        # assume that the area length is twice the area width (so A = L*W = 2W*W = 2W^2)
        buildingWidth = (self.floor.floorArea / 2.) ** 0.5
        manhattanAHUToCornerDistance = 1.5 * buildingWidth  # assume that the AHU is in the middle of the building, and goes from the middle to a corner (so half of width and half of length)
        branchPortion = 0.75  # assume the branch only goes a portion of the way from the AHU to the corner of the building
        branchLineLength = branchPortion * manhattanAHUToCornerDistance
        distanceFromAHUToSensor = 0.6667 * branchLineLength  # static pressure sensor is usually located 2/3rds of the way down the branch line
        pressureDropPerFoot = 0.2 / 100.  # calculate the pressure drop in the duct, per 100 feet per duct
        self.designDPToStaticPressureSensor = distanceFromAHUToSensor * pressureDropPerFoot
        designDischargeStaticPressure = self.designDuctPressure + self.designDPToStaticPressureSensor
        self.designSuctionStaticPressure = -2.5

        self.fanEfficiency = 0.85
        fanBrakeHP = self.designAirflow * (designDischargeStaticPressure - self.designSuctionStaticPressure) / (6356. * self.fanEfficiency)
        flowFactor = 0.67 # corrects for the fricitional factors when moving air, so fan needs to be sized up to correctly handle
        fanSizingSafetyFactor = 1.5

        if fanBrakeHP > 50:
            hpStep = 10.
        elif fanBrakeHP > 10:
            hpStep = 5.
        else:
            hpStep = 1.

        self.fanHP = roundUpToNearest(fanBrakeHP * fanSizingSafetyFactor / flowFactor, hpStep)
        print(self.fanHP, designAirflow)

    def calculateFanTemp(self, timestephours = 0.05):
        #self.calculateEnergyUsage()
        self.electricalPower = ((self.fanSpeed.getActualSpeed() / 100.) ** 3) * self.fanHP * 2544.
        self.fanHeatOutput = (1. - self.fanEfficiency) * self.electricalPower
        if self.dischargeAirflow.getActualAirflow() != 0.:
            self.fanTemp = self.fanHeatOutput / (1.08 * self.dischargeAirflow.getActualAirflow())
        else:
            self.fanTemp = 0.

    def calculateTotalAirflow(self):
        # calculate the total airflow based on the airflow going through child boxes
        totalBoxAirflow = 0.
        for currBox in self.terminalUnitList:
            totalBoxAirflow += currBox.damperAssembly.getActualAirflow()
        self.dischargeAirflow.setActualAirflow(totalBoxAirflow)

    def calculateDesignAirflow(self):
        # go through each of the boxes to determine the design airflow
        totalBoxDesignAirflow = 0.
        for currBox in self.terminalUnitList:
            totalBoxDesignAirflow += currBox.damperAssembly.designAirflow
        return totalBoxDesignAirflow

    def calculateCoolingHeatingRequirements(self):
        # go through each of the rooms to calculate the heating and cooling requirements
        totalDesignCooling = 0
        totalDesignHeating = 0
        for currRoom in self.roomList:
            totalDesignCooling += currRoom.designCoolingLoads
            totalDesignHeating += currRoom.designHeatingLoads
        print("{0} loads: cooling: {1}, heating: {2}".format(self.getName(), round(totalDesignCooling, 0), round(totalDesignHeating, 0)))
        return (totalDesignCooling, totalDesignHeating)

    def calculateFanSpeedEffects(self):
        airflowRatioEffect = (self.dischargeAirflow.getActualAirflow() / self.designAirflow) ** 0.5  # may need to be squared
        speedRatioEffect = (self.fanSpeed.getActualSpeed() / 100.) ** 2.0 #Changing from sqrt to **2
        if airflowRatioEffect < 0.05:
            self.AHU_Duct_Static_Pressure.setActualDP(speedRatioEffect * self.designDPToStaticPressureSensor)
        elif speedRatioEffect < 0.01: # fan speed is less than 10%
            self.AHU_Duct_Static_Pressure.setActualDP(0.)
        else:
            self.AHU_Duct_Static_Pressure.setActualDP(self.designDuctPressure * speedRatioEffect / airflowRatioEffect)
        self.Duct_Static_Pressure_BeforePIU = self.designDuctPressureBeforePIU * (
                    self.AHU_Duct_Static_Pressure.getActualDP() / self.designDuctPressure)
        self.Duct_Static_Pressure_AfterPIU = self.designDuctPressureAfterPIU * (
                    self.AHU_Duct_Static_Pressure.getActualDP() / self.designDuctPressure)

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)

    def createIssueSummary(self):
        return "{" + "'Device': '{0}', 'Is_Broken': {1}, 'Issue_List': {2}".format(self.getName(),
                                                                                   len(self.issueList) > 0,
                                                                                   self.issueList) + "}"

    def resetSimulation(self):
        self.AHU_DAT.setActualTemp(55.)
        self.AHU_RAT.setActualTemp(74.)
        self.AHU_Duct_Static_Pressure.setActualDP(0.2)
        self.fanSpeed.setCommandValue(50.)
        self.CWV.setCommandValue(50.)

    def setOccupancy(self, needHeating=True, needCooling=True):
        self.needHeating = (needHeating != False)
        self.needCooling = (needCooling != False)
        needEither = (needCooling == True or needHeating == True)
        if needEither != self.isOccupied:
            self.isOccupied = needEither
            if self.isOccupied:
                # set the unit to be in occupied mode
                if self.needHeating == True and self.needCooling == False:
                    self.fanSpeed.setCommandValue(20.)
                    self.CWV.setCommandValue(0.)
                else:
                    self.fanSpeed.setCommandValue(20.)
                    self.CWV.setCommandValue(20.)
            else:
                # set the unit to be unoccupied
                self.fanSpeed.setCommandValue(0.)
                self.CWV.setCommandValue(0.)
                self.dischargeAirflow.setActualAirflow(0.)
                self.AHU_Duct_Static_Pressure.setActualDP(0.)
                self.Duct_Static_Pressure_BeforePIU = 0.
                self.Duct_Static_Pressure_AfterPIU = 0.
                self.coolingLoad = 0.

    def logAHUturnOnCooling(self):
        self.numberTerminalUnitsNeedCooling += 1

    def logAHUturnOnHeating(self):
        self.numberTerminalUnitsNeedHeating += 1

    def getOutsideTempAndHumidity(self, outsideConditions):
        self.outdoorTemp = outsideConditions.outdoorTemperature
        self.outdoorHumidity = outsideConditions.outdoorHumidity
        self.OutsideTemp.setActualTemp(self.outdoorTemp)
        self.OutsideHumidity.setActualTemp(self.outdoorHumidity)

    def simulateAHU(self, timeStepHours=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        # update the sensor values
        # set the terminal unit occupancy flag to 0 to count for next iteration
        self.setOccupancy(self.numberTerminalUnitsNeedHeating > 0, self.numberTerminalUnitsNeedCooling > 0)

        self.CHWFlow = self.CWV.getCurrentFlow()
        #self.CHWFlow = self.GPM * (self.CWV.getActualPosition() / 100.)
        if self.isOccupied:
            self.calculateTotalAirflow()  # calculates the airflow through the boxes, which is the airflow through the unit
            self.calculateFanSpeedEffects()  # calculates the duct static pressure reading
            # run any commands on the CHWV, fan speed, dampers, etc, and then update the sensor values
            self.CWRT = min(self.CWstartTemp + self.designDTacrosscoil, self.AHU_RAT.getActualTemp() - self.RATdeltaT)
            self.coolingLoad = 500. * (self.CWRT - self.CWstartTemp) * self.CHWFlow
            DAT = self.AHU_RAT.getActualTemp() - (
                        self.coolingLoad / (max(1., self.dischargeAirflow.getActualAirflow()) * 1.08))
            if DAT < self.CWstartTemp + self.CHWstartDeltaT:
                self.AHU_DAT.setActualTemp(self.CWstartTemp + self.CHWstartDeltaT)
                self.coolingLoad = self.dischargeAirflow.getActualAirflow() * 1.08 * (
                            self.AHU_RAT.getActualTemp() - (self.CWstartTemp + self.CHWstartDeltaT))
                self.CWRT = self.CWstartTemp + (self.coolingLoad / (max(1., self.CHWFlow) * 500.))
            else:
                self.AHU_DAT.setActualTemp(DAT)

            deltaCWVpos = (self.coolingOutput.calculateCWVCoolingOutput(self.AHU_DAT.getSensorReading() - self.DATSP))
            if self.needHeating == True and self.needCooling == False:
                self.CWV.setCommandValue(0.)
            else:
                self.CWV.setCommandValue(self.CWV.getSensorReading() + deltaCWVpos)
            # command fan speed to keep duct static pressure at setpoint
            deltaFanSpeed = self.fanspeedoutput.calculatefanspeedOutput(self.AHU_Duct_Static_Pressure.getSensorDP() - self.StaticPressureSP)
            self.fanSpeed.setCommandValue(self.fanSpeed.getSensorReading() - deltaFanSpeed)
        else:
            self.CWRT = self.CWstartTemp
            self.AHU_DAT.setActualTemp(self.AHU_RAT.getActualTemp())
        self.numberTerminalUnitsNeedHeating = 0
        self.numberTerminalUnitsNeedCooling = 0

    def addTerminalUnit(self, newTerminalUnit):
        self.terminalUnitList.append(newTerminalUnit)

    def addRoom(self, newRoom):
        self.roomList.append(newRoom)

    def addPoint(self, pointInfo):
        pointName = pointInfo.getName()
        if not (pointName in self.pointNameList):
            self.pointsList.append(pointInfo)
            self.pointNameList.append(pointName)

    def randomizeSensorStatus(self):
        for point in self.pointsList:
            point.randomizeSensorStatus()

    def getName(self):
        return self.AHU_Name

    def getEquipmentType(self):
        return "AHU"

    def setUpOutsideTemperature(self, hasSensor=True, setUpWorking=None):
        self.OutsideTemp = temperatureSensor(self, "Outside Temperature", logValue=(hasSensor == True), setUpWorking=setUpWorking)
        self.OutsideTemp.probabilityOfSensor = 0.7
        self.addPoint(self.OutsideTemp)

    def setUpOutsideHumidity(self, hasSensor=True, setUpWorking=None):
        self.OutsideHumidity = temperatureSensor(self, "Outside Humidity", logValue=(hasSensor == True), setUpWorking=setUpWorking)
        self.OutsideHumidity.probabilityOfSensor = 0.7
        self.addPoint(self.OutsideHumidity)

    def setUpCWVOutput(self, hasSensor=True, setUpWorking=None):
        self.CWV = CWVOutput(parentDevice=self, outputName="CWV Pos", logValue=(hasSensor == True), setUpWorking=setUpWorking)
        self.CWV.probabilityOfSensor = 0.7
        self.addPoint(self.CWV)

    def setUpDATSensor(self, hasSensor=True, setUpWorking=None):
        self.AHU_DAT = temperatureSensor(self, "AHU DAT Temp", logValue=(hasSensor == True), setUpWorking=setUpWorking)
        self.AHU_DAT.probabilityOfSensor = 0.7
        self.addPoint(self.AHU_DAT)

    def setUpAHUDATSP(self, hasSensor=True, setUpWorking=None):
        self.AHUDATSP = AHUDATSP(parentDevice=self, sensorName="AHU DAT SP", logValue=(hasSensor==True),
                                               setUpWorking=setUpWorking)
        self.AHUDATSP.probabilityOfSensor = 0.7
        self.addPoint(self.AHUDATSP)
        self.AHUDATSP.setAHUDATSP(self.DATSP)

    def setUpRATSensor(self, hasSensor=True, setUpWorking=None):
        self.AHU_RAT = temperatureSensor(self, "AHU RAT Temp", logValue=(hasSensor == True), setUpWorking=setUpWorking)
        self.AHU_RAT.probabilityOfSensor = 0.7
        self.addPoint(self.AHU_RAT)

    def setUpStaticPressureSensor(self, hasSensor=True, setUpWorking=None):
        self.AHU_Duct_Static_Pressure = airPressureTransducer(self, "AHU SAP", designDP=0.2,
                                                              logValue=(hasSensor == True), setUpWorking=setUpWorking)
        self.AHU_Duct_Static_Pressure.setPressureRange(minReading=0., maxReading=4.)
        self.AHU_Duct_Static_Pressure.probabilityOfSensor = 0.7
        self.addPoint(self.AHU_Duct_Static_Pressure)

    def setUpStaticPressureSPsensor(self, hasSensor=True, setUpWorking=None):
        self.SPSP = StaticPressureSP(parentDevice=self, sensorName="Static Pressure SP", logValue=(hasSensor==True),
                                               setUpWorking=setUpWorking)
        self.SPSP.probabilityOfSensor = 0.7
        self.addPoint(self.SPSP)
        self.SPSP.setStaticPressureSP(self.StaticPressureSP)

    def setUpDischargeAirflowSensor(self, hasSensor=True, setUpWorking=None):
        self.dischargeAirflow = airflowStation(self, "Supply Airflow", logValue=(hasSensor == True),
                                               setUpWorking=setUpWorking)
        self.dischargeAirflow.probabilityOfSensor = 0.7
        self.addPoint(self.dischargeAirflow)

    def setUpFanSpeed(self, isVariableSpeed=True, hasSensor=True, setUpWorking=None):
        self.fanSpeed = speedControl(self, "Fan Speed", isVariableSpeed=isVariableSpeed, logValue=(hasSensor == True),
                                     setUpWorking=setUpWorking)
        self.fanSpeed.probabilityOfSensor = 0.7
        self.addPoint(self.fanSpeed)

    def logConditions(self, timestamp, logAsJSON=True):
        isJSONFormat = (logAsJSON == True)
        return logDevice(timestamp, self, isJSONFormat, logIssueStatus=True)

    def calculateFloorEnergyUsage(self, timestephours=0.05):
        roomEnergy = 0.
        for room in self.roomList:
            roomEnergy += room.calculateEnergyUsage()
        fanEnergy = (self.fanHP * (self.fanSpeed.getActualSpeed() / 100.) ** 3.) * 745.7   # in Watts
        return fanEnergy + roomEnergy

class ChilledWaterSystem:
    def __init__(self, parentBuilding):
        self.parentBuilding = parentBuilding
        self.PumpSystem = PumpSystem(parentBuilding=parentBuilding)
        self.pointsList = []
        self.issueList = []
        self.AirHandlerList = []
        self.CWstartTemp = 45.
        self.RATdeltaT = 3.
        self.CHWstartDeltaT = 3.
        self.designDTacrosscoil = 10.

    def getName(self):
        return "Chilled Water System"

    def getEquipmentType(self):
        return "CHW System"

    def addAirHandler(self, newAirHandler):
        self.AirHandlerList.append(newAirHandler)
        self.PumpSystem.addAirHandler(newAirHandler)

    def calculate(self):
        currentGPM = 0.
        totalTempWeighted = 0.
        for AHU in self.AirHandlerList:
            currentGPM += AHU.CHWFlow
            totalTempWeighted += AHU.CHWFlow * AHU.CWRT
        if currentGPM > 0.:
            CHWreturnTemp = totalTempWeighted / currentGPM
        else:
            CHWreturnTemp = self.CWstartTemp

        self.PumpSystem.calculateDP()
        return CHWreturnTemp

    def createIssueSummary(self):
        summaryList = ""
        summaryList = addLogToRunningList(summaryList, "{" + "'Device': '{0}', 'Is_Broken': {1}, 'Issue_List': {2}".format(self.getName(),
                                                                                   len(self.issueList) > 0,
                                                                                   self.issueList) + "}", True)
        summaryList = addLogToRunningList(summaryList, self.PumpSystem.createIssueSummary(), True)
        return summaryList


    def sizePlant(self):
        self.PumpSystem.setDesignDP(14.)
        self.PumpSystem.sizePumpSystem()
        self.designFlow = self.PumpSystem.designFlow
        self.designTonnage = self.designFlow * (self.designDTacrosscoil)


    def logConditions(self, timestamp, logAsJSON=True):
        isJSONFormat = (logAsJSON == True)
        logList = ""
        logList = addLogToRunningList(logList, logDevice(timestamp, self, isJSONFormat, logIssueStatus=True), isJSONFormat)
        logList = addLogToRunningList(logList, self.PumpSystem.logConditions(timestamp, isJSONFormat))
        return logList

class PumpSystem:
    def __init__(self, parentBuilding):
        self.parentBuilding = parentBuilding
        self.pointsList = []
        self.pointNameList = []
        self.CHWVList = []
        self.PumpList = []
        self.addPumps(setUpWorking=True)
        self.setUpDPSensor(setUpWorking=True)
        self.PumpOutputLoop = coolingPumpOutput()
        self.setDesignDP(14.)
        self.setDPSetpoint(self.designDP)
        self.designFlow = 100.
        self.oldDP = self.designDP
        self.oldGPM = 30.
        self.oldSpeed = 100.

        self.issueList = []

    def sizePumpSystem(self):
        designFlow = 0
        for CWV in self.CHWVList:
            designFlow += CWV.designGPM
        self.designFlow = designFlow
        self.oldGPM = self.designFlow
        self.oldSpeed = 100.
        return self.designFlow

    def addAirHandler(self, newAirHandler):
        self.CHWVList.append(newAirHandler.CWV)

    def addPumps(self, hasSensor=True, setUpWorking=None):
        self.pump = PumpOutput(parentDevice=self, outputName="Pump Pos", logValue=(hasSensor == True), setUpWorking=setUpWorking)
        self.pump.setCommandValue(50.)
        self.addPoint(self.pump)

    def setUpDPSensor(self, hasSensor=True, setUpWorking=True):
        self.systemDP = waterSystemDP(parentDevice=self, sensorName="CHWS DP", logValue=hasSensor, setUpWorking=setUpWorking)
        self.addPoint(self.systemDP)
        self.dpSetpoint = waterSystemDP(parentDevice=self, sensorName="DP Stpt", logValue=hasSensor,
                                      setUpWorking=True)
        self.dpSetpoint.analogInputPoint.readingUncertainty = 0. # make the setpoint have no uncertainty or offset
        self.dpSetpoint.analogInputPoint.readingOffset = 0.
        self.addPoint(self.dpSetpoint)

    def addPoint(self, pointInfo):
        pointName = pointInfo.getName()
        if not (pointName in self.pointNameList):
            self.pointsList.append(pointInfo)
            self.pointNameList.append(pointName)

    def setDesignDP(self, newDesignDP):
        self.systemDP.setActualDP(newDesignDP) # so that is changed to a non-zero value, and tracks the system
        self.designDP = newDesignDP
        for CWV in self.CHWVList:
            CWV.setDesignDP(self.designDP)

    def setDPSetpoint(self, dpSetpoint):
        self.systemDP.setActualDP(dpSetpoint)  # so that is changed to a non-zero value, and tracks the system
        self.dpSetpoint.setActualDP(dpSetpoint)

    def getName(self):
        return "Pump List"

    def getEquipmentType(self):
        return "Pump"

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)

    def createIssueSummary(self):
        return "{" + "'Device': '{0}', 'Is_Broken': {1}, 'Issue_List': {2}".format(self.getName(),
                                                                                   len(self.issueList) > 0,
                                                                                   self.issueList) + "}"

    def calculateDP(self):
        currentGPM = 0.
        actualDP = self.systemDP.getActualDP()
        readingDP = self.systemDP.getSensorDP()
        DPSP = self.dpSetpoint.getActualDP()
        for CWV in self.CHWVList:
            CWV.setCurrentDP(actualDP)
            #AHU.setCurrentDP(14.)
            currentGPM += CWV.getCurrentFlow()
        currentGPM = max(0.1, currentGPM)
        gpmRatio = self.oldGPM / max(0.1, currentGPM)
        self.currentDP = min(50., (gpmRatio ** 2) * self.oldDP)
        deltaPumpPos = (self.PumpOutputLoop.calculatePumpCoolingOutput(DPSP - readingDP))
        print("{0}\t{1}\t{2}\t{3}\t{4}".format(deltaPumpPos, self.pump.getActualPosition(), actualDP, currentGPM, self.CHWVList[0].getActualPosition()))
        self.pump.setCommandValue(self.pump.getActualPosition() + deltaPumpPos)
        pumpSpeedEffects = (((self.pump.getActualPosition()) / max(1., self.oldSpeed)) ** 2.)
        print(self.pump.getActualPosition(), self.oldSpeed, pumpSpeedEffects)
        actualDP *= pumpSpeedEffects
        self.systemDP.setActualDP(actualDP)


        currentGPM = 0
        for CWV in self.CHWVList:
            currentGPM += CWV.setCurrentDP(actualDP)
            # AHU.CWV.setFlowFactor(gpmRatio)
            #AHU.setGPMratio(1.)

        self.oldDP = max(0.1, actualDP)
        self.oldGPM = currentGPM
        self.oldSpeed = max(0.1, self.pump.getActualPosition())

    def logConditions(self, timestamp, logAsJSON=True):
        isJSONFormat = (logAsJSON == True)
        return logDevice(timestamp, self, isJSONFormat, logIssueStatus=True)


class BoilerSystem:
    def __init__(self, parentBuilding):
        self.parentBuilding = parentBuilding
        self.TerminalUnitList = []
        self.HWstartTemp = 90.
        self.RATdeltaT = 3.
        self.HWstartDeltaT = 3.
        self.designDTacrosscoil = 20.

    def addTerminalUnit(self, newTerminalUnit):
        self.TerminalUnitList.append(newTerminalUnit)

    def calculate(self):
        totalGPM = 0.
        totalTempWeighted = 0.
        for TU in self.TerminalUnitList:
            totalGPM += TU.hotWaterReturnTemperatureandFlow()[1]
            totalTempWeighted += TU.hotWaterReturnTemperatureandFlow()[1] * TU.hotWaterReturnTemperatureandFlow()[0]
        if totalGPM > 0.:
            HWRT = totalTempWeighted / totalGPM
        else:
            HWRT = self.HWstartTemp
        return HWRT

class RoomRepresentation:
    def __init__(self,
                 roomName,
                 parentBuilding,
                 parentFloor,
                 parentAHU,
                 isRoof=False,
                 randomizeSensors=True,
                 setUpWorking=None):
        self.parentBuilding = parentBuilding
        self.parentFloor = parentFloor
        self.parentAHU = parentAHU
        self.randomizeSensors = randomizeSensors
        self.setUpWorking = setUpWorking
        self.parentAHU.addRoom(self)
        self.isRoof = (isRoof == True)

        self.roomName = roomName
        self.issueList = []
        self.pointsList = []
        self.pointNameList = []

        self.setDimensions()
        self.setOccupants()
        self.setLightsDesign()
        self.setBaseLoads()
        self.setOARequirement()
        self.roomTempStpts = self.parentBuilding.roomTempStpts.copyRoomSetpoints()

        self.occupancyInfo = RoomUsageStatus(self.maxLightHeatLoad, self.maxOccupantHeatLoad, self.baseHeatLoad,
                                             self.parentBuilding.occupancySchedule)  # all of these are in kW
        self.setDesignRequirements()
        self.terminalUnit = None
        self.fitTerminalUnit()

    def fitTerminalUnit(self):
        self.terminalUnit = terminalUnitBox(boxName=self.roomName + " TU",
                                            parentAHU=self.parentAHU,
                                            roomServed=self,
                                            currentFloor=self.parentFloor,
                                            parentBuilding=self.parentBuilding,
                                            randomizeSensors=self.randomizeSensors,
                                            setUpWorking=self.setUpWorking)

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)

    def addPoint(self, pointInfo):
        pointName = pointInfo.getName()
        if not (pointName in self.pointNameList):
            self.pointsList.append(pointInfo)
            self.pointNameList.append(pointName)

    def getName(self):
        return self.roomName

    def getEquipmentType(self):
        return "Room"

    def logConditions(self, timestamp, logAsJSON=True):
        isJSONFormat = (logAsJSON == True)
        # log all of the room information
        logString = logDevice(timestamp, self, isJSONFormat, logIssueStatus=True)

        # log the terminal unit information
        logString = addLogToRunningList(logString, self.terminalUnit.logConditions(timestamp, isJSONFormat),
                                        isJSONFormat)

        return logString

    def setDimensions(self, roomWidth=None, roomLength=None, roomHeight=None, eastOrWestOuterWall=None,
                      northOrSouthOuterWall=None):
        if roomWidth is None:
            roomWidth = round(selectRandomFromRangeList(((5, 15), (15, 30), (30, 80)), (0.1, 0.5, 0.4)), 1)
        if roomLength is None:
            roomLength = round(selectRandomFromRangeList(((5, 15), (15, 30), (30, 80)), (0.1, 0.5, 0.4)), 1)
        if roomHeight is None:
            roomHeight = round(selectRandomElement((9, 11, 14), (0.7, 0.1, 0.2)), 1)  # room height in feet
        if eastOrWestOuterWall is None:
            eastOrWestOuterWall = selectRandomElement(("East", "West", "Inside"), (0.2, 0.2, 0.6))
        if northOrSouthOuterWall is None:
            northOrSouthOuterWall = selectRandomElement(("North", "South", "Inside"), (0.2, 0.2, 0.6))

        self.roomWidth = roomWidth
        self.roomLength = roomLength
        self.roomHeight = roomHeight
        self.floorSize = round(self.roomWidth * self.roomLength, 0)
        self.roomVolume = round(self.roomWidth * self.roomLength * self.roomHeight, 0)  # room volume in cubic feet
        self.roomHeatInertia = self.roomVolume * 0.075 * 0.24  # convert the room volume to pounds (0.075 lbs / ft^3), and then to Btu's to increase 1 degree (0.24 degF / (BTU * lb))
        self.wallFaceTypes = "E/W: {0}, N/S: {1}, Roof: {2}".format(eastOrWestOuterWall, northOrSouthOuterWall, self.isRoof)

        # create the walls, and roof
        self.roomFaceList = []
        self.roomFaceList.append(RoomFace(self.roomWidth * self.roomHeight,
                                          self.parentBuilding.constructionMaterialParameters,
                                          degreeAngleFromHorizontal=0,
                                          degreeAngleAzimuth=90,
                                          isExterior=(eastOrWestOuterWall == "East"),
                                          faceName="East Wall"))
        self.roomFaceList.append(RoomFace(self.roomWidth * self.roomHeight,
                                          self.parentBuilding.constructionMaterialParameters,
                                          degreeAngleFromHorizontal=0,
                                          degreeAngleAzimuth=270,
                                          isExterior=(eastOrWestOuterWall == "West"),
                                          faceName="West Wall"))
        self.roomFaceList.append(RoomFace(self.roomLength * self.roomHeight,
                                          self.parentBuilding.constructionMaterialParameters,
                                          degreeAngleFromHorizontal=0,
                                          degreeAngleAzimuth=0,
                                          isExterior=(northOrSouthOuterWall == "North"),
                                          faceName="North Wall"))
        self.roomFaceList.append(RoomFace(self.roomLength * self.roomHeight,
                                          self.parentBuilding.constructionMaterialParameters,
                                          degreeAngleFromHorizontal=0,
                                          degreeAngleAzimuth=180,
                                          isExterior=(northOrSouthOuterWall == "South"),
                                          faceName="South Wall"))
        self.roomFaceList.append(RoomFace(self.roomLength * self.roomWidth,
                                          self.parentBuilding.constructionMaterialParameters,
                                          degreeAngleFromHorizontal=90,
                                          degreeAngleAzimuth=0,
                                          isExterior=self.isRoof,  # default to whatever the floor value is
                                          faceName="Ceiling"))
        self.roomFaceList.append(RoomFace(self.roomLength * self.roomWidth,
                                          self.parentBuilding.constructionMaterialParameters,
                                          degreeAngleFromHorizontal=-90,
                                          degreeAngleAzimuth=0,
                                          isExterior=False,
                                          faceName="Floor"))

    def setBaseLoads(self):
        self.baseHeatLoad = round(selectRandomFromRangeList(((0, 1), (1, 2), (2, 5)), (0.5, 0.3, 0.2)),
                                  2)  # base heat load in kW
        self.baseLoadSafetyFactor = round(
            selectRandomFromRangeList(((0.5, 0.9), (0.9, 1.1), (1.1, 1.3), (1.3, 2)), (0.05, 0.05, 0.7, 0.2)), 2)
        self.designBaseHeatLoad = self.baseHeatLoad * self.baseLoadSafetyFactor

    def setLightsDesign(self):
        self.lightPerArea = round(
            selectRandomFromRangeList(((0.5, 0.9), (0.9, 1.1), (1.1, 1.3), (1.3, 2)), (0.4, 0.05, 0.5, 0.05)),
            2)  # light in W/sqft
        self.lightPercentInSpace = 0.8  # percent of light heat that goes into space
        self.plenumLightLoadFactor = (
                                                 1. - self.lightPercentInSpace) / self.lightPercentInSpace  # given the light in the space, use this factor to calculate the amount of light going to the plenum
        self.lightingSafetyFactor = round(
            selectRandomFromRangeList(((0.5, 0.9), (0.9, 1.1), (1.1, 1.3), (1.3, 2)), (0.05, 0.05, 0.7, 0.2)), 2)
        self.designLightingPower = self.floorSize * self.lightPerArea
        self.maxLightHeatLoad = round(self.lightPercentInSpace * self.floorSize * self.lightPerArea / 1000,
                                      2)  # lighting heat load in kW
        self.designLightHeatLoad = self.maxLightHeatLoad * self.lightingSafetyFactor  # lighting heat load in kW, with safety factor

    def setOccupants(self):
        self.occupantDensity = round(
            selectRandomFromRangeList(((0, 0), (0, 4), (4, 6), (6, 15), (15, 30)), (0.05, 0.1, 0.6, 0.15, 0.1)),
            2) / 1000.  # number of occupants / 1000 sqft
        self.numOccupants = round(self.occupantDensity * self.floorSize, 0)
        self.heatFromOccupant = 100.  # heat in watts
        self.maxOccupantHeatLoad = self.numOccupants * self.heatFromOccupant / 1000.  # occupant heat load in kW
        self.occupantSafetyFactor = round(
            selectRandomFromRangeList(((0.5, 0.9), (0.9, 1.1), (1.1, 1.3), (1.3, 2)), (0.05, 0.05, 0.7, 0.2)), 2)
        self.designOccupantHeatLoad = self.occupantSafetyFactor * self.maxOccupantHeatLoad

    def setOARequirement(self):
        self.OAperUnitArea = selectRandomElement((0, 0.06, 0.12, 0.18), (0.05, 0.6, 0.25, 0.1))  # cfm OA per sq ft
        self.OAperPerson = 5.  # cfm OA per occupant
        self.minOA = self.OAperUnitArea * self.floorSize + self.OAperPerson * self.numOccupants  # minimum CFM for the room

    def setDesignRequirements(self):
        # gets the required cooling (in BTU/h) required per the design
        internalCoolingLoads = (
                                           self.designBaseHeatLoad + self.designLightHeatLoad + self.designOccupantHeatLoad) * 3412  # in kW, so convert to Btu/hr
        externalCoolingLoads = 0
        for currWall in self.roomFaceList:
            externalCoolingLoads += currWall.getHeatFlowToRoom(self.roomTempStpts.designOccCoolingStpt, self.parentFloor.averageAreaTemp,
                                                               self.parentBuilding.coolingOutsideParameters)  # returns in Btu/hr
        externalCoolingLoads *= 1.2  # add a safety factor
        self.designCoolingLoads = internalCoolingLoads + externalCoolingLoads

        # gets the required heating (in BTU/h) required per the design
        internalHeatingLoads = 0  # assume nothing is on
        externalHeatingLoads = 0
        for currWall in self.roomFaceList:
            externalHeatingLoads += currWall.getHeatFlowToRoom(self.roomTempStpts.designOccHeatingStpt, self.parentFloor.averageAreaTemp,
                                                               self.parentBuilding.heatingOutsideParameters)  # returns in Btu/hr
        externalHeatingLoads *= 1.2  # add a safety factor
        self.designHeatingLoads = -(internalHeatingLoads + externalHeatingLoads)  # negate to give a positive value

    def resetSimulation(self):
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        self.terminalUnit.resetSimulation()

    def resetSimulationTime(self, hourOfDay=0.):
        # reset all references to time within the simulation to a new time
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        self.occupancyInfo.resetToTime(hourOfDay)
        self.terminalUnit.setOccupancy(self.occupancyInfo.isOccupied)
        self.terminalUnit.resetSimulationTime(hourOfDay)

    def getHeatFlowThroughFaces(self, outdoorConditions):
        # outputs in Btu / hr
        faceHeatFlow = 0
        for currWall in self.roomFaceList:
            # the face determines whether to use the average floor temp or the outside temp based on if the wall is exterior
            faceHeatFlow += currWall.getHeatFlowToRoom(self.terminalUnit.roomTemp.getActualTemp(), self.parentFloor.averageAreaTemp,
                                                       outdoorConditions)  # returns in Btu/hr
        return faceHeatFlow

    def getHeatFlowIntoRoom(self):
        # get the net heat flow (in BTU/hr) INTO the room, and the net airflow INTO the room
        # get the impact of the terminal unit
        unitAirflow, unitAirTemp = self.terminalUnit.runAirflowSimulation()
        terminalUnitLoad = (unitAirTemp - self.terminalUnit.roomTemp.getActualTemp()) * 1.08 * unitAirflow
        # get the impact of the internal loads
        internalLoads = self.occupancyInfo.calculateLoads() * 3412.  # occupancy internal loads are calculated in kW, so convert to Btu/hr

        # get the impact of the external loads
        externalLoads = self.getHeatFlowThroughFaces(
            self.parentBuilding.outsideConditions)  # external loads returned in Btu/hr
        return unitAirflow, terminalUnitLoad + internalLoads + externalLoads

    def getHeatFlowIntoPlenum(self):
        return self.plenumLightLoadFactor * self.occupancyInfo.getLightLoad() * 3412.  # convert from kW to BTU/hr

    def simulateRoom(self, timeStepHours=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        # move both the occupancy and the terminal unit forward by the time step
        self.occupancyInfo.incrementTime(timeStepHours)
        self.terminalUnit.setOccupancy(self.occupancyInfo.isOccupied)
        self.terminalUnit.runControlSequence(timeStepHours, self.occupancyInfo.getAHUstatus())

        airFlowIntoRoom, heatFlowIntoRoom = self.getHeatFlowIntoRoom()
        self.terminalUnit.roomTemp.setActualTemp(
            heatFlowIntoRoom * timeStepHours / self.roomHeatInertia + self.terminalUnit.roomTemp.getActualTemp())
        self.calculateEnergyUsage()

        return airFlowIntoRoom, self.terminalUnit.roomTemp.getActualTemp()

    def calculateEnergyUsage(self):
        baseLoad = self.occupancyInfo.getBaseLoad()  # assumed that the base load is coming from electronics or other power supply
        lightLoad = self.designLightingPower * self.occupancyInfo.getLightingFactor()
        TUFanandHeaterLoad = self.terminalUnit.calculateEnergyUsage()
        return baseLoad + lightLoad + TUFanandHeaterLoad

class terminalUnitBox:
    # class which describes a terminal unit
    def __init__(self, boxName="New Terminal Unit",
                 parentAHU=None,
                 roomServed=None,
                 currentFloor=None,
                 parentBuilding=None,
                 randomizeSensors=True,
                 setUpWorking=None):
        self.issueList = []
        self.pointsList = []
        self.pointNameList = []

        # set all references to the AHU
        self.setAHU(parentAHU)

        # set all references to the room served
        self.setRoom(roomServed)

        # set all references to the floor and building
        self.setFloor(currentFloor)
        self.setBuilding(parentBuilding)
        self.assignHeaterFanPackageType()
        self.boxName = boxName + '_' + self.unitType

        # set the control sequence function
        self.coolingOutput = coolingAirflowOutput()
        self.heatingOutput = heatingOutput()
        self.heatingValveOutput = coolingCWVOutput()

        # initialize the occupancy
        self.setUpOccupancyStatus(setUpWorking=setUpWorking)
        self.setOccupancy(occupancyStatus=True)
        self.setDesignAirflows()

        self.setUpAirsideSystem(pressureDrop=self.AHU.designDuctPressureBeforePIU - self.AHU.designDuctPressureAfterPIU,
                                setUpWorking=setUpWorking)
        self.setUpNeedCoolingStatus(setUpWorking=setUpWorking)
        self.setUpNeedHeatingStatus(setUpWorking=setUpWorking)
        self.setNeedHeat()
        self.setNeedCool()

        # set up all of the applicable actuators and sensors
        self.setUpRoomCoolingStpt(setUpWorking=setUpWorking)
        self.setUpRoomHeatingStpt(setUpWorking=setUpWorking)
        self.setUpDATSensor(hasSensor=True, setUpWorking=setUpWorking)
        self.setUpRoomTemp(hasSensor=True, setUpWorking=setUpWorking)
        # self.assignHeaterFanPackageType(TUHeatingType=None)
        self.setUpHeaterFanPackage(setUpWorking=setUpWorking)

        if randomizeSensors:
            self.randomizeSensorStatus()
        self.dTacrosscoil = 20.

    def setAHU(self, AHU=None):
        self.AHU = AHU
        self.addPoint(self.AHU.AHU_DAT)
        self.AHU.addTerminalUnit(self)

    def setBuilding(self, parentBuilding = None):
        self.parentBuilding = parentBuilding
        self.applicableUnitTypes = parentBuilding.defaultTUHeatingType

    def setRoom(self, roomServed=None):
        self.roomServed = roomServed
        self.roomStpts = self.roomServed.roomTempStpts

    def setFloor(self, currentFloor=None):
        self.currentFloor = currentFloor
        # terminal unit should already be associated when the floor is set up

    def setDesignAirflows(self):
        designCoolingBTU = self.roomServed.designCoolingLoads
        designHeatingBTU = self.roomServed.designHeatingLoads
        heatingDAT = self.AHU.design_TU_Heating_DAT
        coolingDAT = self.AHU.design_AHU_Cooling_DAT
        coolingStpt = self.roomStpts.designOccCoolingStpt
        heatingStpt = self.roomStpts.designUnoccHeatingStpt  # peak heating load should be in the middle of the night, when unoccupied
        self.maxCoolingCFM = roundUpToNearest(designCoolingBTU / (1.08 * (coolingStpt - coolingDAT)), 10.)
        self.minCoolingCFM = 0. # roundToNearest(self.maxCoolingCFM * 0.1, 10)
        self.maxHeatingCFM = roundUpToNearest(max(100., designHeatingBTU / (1.08 * (heatingDAT - heatingStpt))), 10)
        additionalHeatingRequired = self.maxHeatingCFM * (heatingStpt - coolingDAT) * 1.08
        self.roomServed.designHeatingLoads += additionalHeatingRequired # add in the heat required to get the air coming from the AHU up to room temperature (assuming that the AHU is in cooling mode)
        self.minHeatingCFM = roundUpToNearest(max(50, self.maxHeatingCFM * 0.1), 10)
        self.coolingOutput.setMinMaxAirflow(self.minCoolingCFM, self.maxCoolingCFM)

    def assignHeaterFanPackageType(self, TUHeatingType=None):
        TUHeatingTypeList = ["VAV", "VAVR", "VAVRE", "PIU", "SPIU"]
        removeList = list(self.applicableUnitTypes)
        if TUHeatingType not in TUHeatingTypeList:
            TUHeatingType = None
        heaterCutoff = 0.5 * 3412  # BTU's per 1 kW
        shouldHaveHeater = self.roomServed.designHeatingLoads >= heaterCutoff
        hasHeater = selectRandomElement([shouldHaveHeater, not shouldHaveHeater], [0.9, 0.1])
        if hasHeater is False:
            TUHeatingType=self.applicableUnitTypes[1]
        else:
            TUHeatingType = selectRandomElement([self.applicableUnitTypes[0],selectRandomElement(list(set(TUHeatingTypeList).difference(set(removeList))))],[0.9,0.1])
        self.unitType = TUHeatingType
        if self.unitType == "VAVR":
            self.parentBuilding.BoilerSystem.addTerminalUnit(self)

    def setOccupancy(self, occupancyStatus=True):
        occupancyStatus = (occupancyStatus != False) + 0
        if occupancyStatus != self.occupancyStatus.getActualOccupancyStatus():
            self.occupancyStatus.setActualOccupancyStatus(occupancyStatus)
            if self.occupancyStatus.getActualOccupancyStatus() == 1:
                self.heatingStpt = self.roomStpts.occupiedHeatingStpt
                self.coolingStpt = self.roomStpts.occupiedCoolingStpt
            else:
                self.heatingStpt = self.roomStpts.unoccupiedHeatingStpt
                self.coolingStpt = self.roomStpts.unoccupiedCoolingStpt

    def setNeedHeat(self):
        if self.needHeatingStatus.getActualNeedHeatingStatus() == 1:
            self.AHU.logAHUturnOnHeating()

    def setNeedCool(self):
        if self.needCoolingStatus.getActualNeedCoolingStatus() == 1:
            self.AHU.logAHUturnOnCooling()


    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)

    def createIssueSummary(self):
        return "{" + "'Device': '{0}', 'Is_Broken': {1}, 'Issue_List': {2}".format(self.boxName,
                                                                                   len(self.issueList) > 0,
                                                                                   self.issueList) + "}"

    def addPoint(self, pointInfo):
        pointName = pointInfo.getName()
        if not (pointName in self.pointNameList):
            self.pointsList.append(pointInfo)
            self.pointNameList.append(pointName)

    def randomizeSensorStatus(self):
        for point in self.pointsList:
            point.randomizeSensorStatus()

    def getName(self):
        return self.boxName

    def getEquipmentType(self):
        return self.unitType

    def logConditions(self, timestamp, logAsJSON=True):
        isJSONFormat = (logAsJSON == True)
        return logDevice(timestamp, self, isJSONFormat, logIssueStatus=True)

    def resetSimulation(self):
        self.damperAssembly.reset()
        self.roomTemp.setActualTemp(76.)  # initialize at 70 degF

    def resetSimulationTime(self, hourOfDay=0.):
        # reset all references to time within the simulation to a new time
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        pass  # no time records within the room

    def setUpHeaterFanPackage(self, setUpWorking=None):
        self.heaterFanPackage = heaterFanPackage(self.maxHeatingCFM*0.75, self.roomServed.designHeatingLoads, self.unitType,
                                                 parentDevice=self, outputName="TU Temp",
                                                setUpWorking=setUpWorking)
        self.heaterFanPackage.probabilityOfSensor = 0.7
        self.addPoint(self.heaterFanPackage)

    def setUpAirsideSystem(self, pressureDrop=0.3, setUpWorking=None):
        maxBoxCFM = max(self.maxCoolingCFM, self.maxHeatingCFM)
        self.damperAssembly = damperAirflow(maxBoxCFM, pressureDrop, parentDevice=self, outputName="Airflow",
                                            isVariableVolume=True, logValue=True, setUpWorking=setUpWorking)
        self.damperAssembly.probabilityOfSensor = 0.7
        self.damperPos = self.damperAssembly.damperPosition
        self.addPoint(self.damperAssembly)

    def setUpDATSensor(self, hasSensor=True, setUpWorking=None):
        self.dischargeTemp = temperatureSensor(parentDevice=self, sensorName="TU DAT", logValue=hasSensor,
                                               setUpWorking=setUpWorking)
        self.dischargeTemp.probabilityOfSensor = 0.6
        self.addPoint(self.dischargeTemp)

    def setUpRoomTemp(self, hasSensor=True, setUpWorking=None):
        self.roomTemp = temperatureSensor(parentDevice=self, sensorName="Room Temp", logValue=hasSensor,
                                               setUpWorking=setUpWorking)
        self.roomTemp.probabilityOfSensor = 0.95
        self.addPoint(self.roomTemp)

    def setUpOccupancyStatus(self, hasSensor=True, setUpWorking=None):
        self.occupancyStatus = OccupancyStatus(parentDevice=self, sensorName="Occupancy Status", logValue=hasSensor,
                                               setUpWorking=setUpWorking)
        self.occupancyStatus.probabilityOfSensor = 0.7
        self.addPoint(self.occupancyStatus)

    def setUpNeedCoolingStatus(self, hasSensor=False, setUpWorking=None):
        self.needCoolingStatus = NeedCoolingStatus(parentDevice=self, sensorName="Need Cooling Status",
                                                                     logValue=hasSensor, setUpWorking=setUpWorking)
        self.needCoolingStatus.probabilityOfSensor = 0.7
        self.addPoint(self.needCoolingStatus)

    def setUpNeedHeatingStatus(self, hasSensor=False, setUpWorking=None):
        self.needHeatingStatus = NeedHeatingStatus(parentDevice=self, sensorName="Need Heating Status",
                                                                     logValue=hasSensor, setUpWorking=setUpWorking)
        self.needHeatingStatus.probabilityOfSensor = 0.7
        self.addPoint(self.needHeatingStatus)

    def setUpRoomCoolingStpt(self, hasSensor=True, setUpWorking=None):
        self.roomCoolingStpt = RoomCoolingStpt(parentDevice=self, sensorName="Room Cooling Stpt", logValue=hasSensor,
                                               setUpWorking=setUpWorking)
        self.roomCoolingStpt.probabilityOfSensor = 0.7
        self.addPoint(self.roomCoolingStpt)

    def setUpRoomHeatingStpt(self, hasSensor=True, setUpWorking=None):
        self.roomHeatingStpt = RoomCoolingStpt(parentDevice=self, sensorName="Room Heating Stpt", logValue=hasSensor,
                                               setUpWorking=setUpWorking)
        self.roomHeatingStpt.probabilityOfSensor = 0.7
        self.addPoint(self.roomHeatingStpt)

    def runControlSequence(self, timeStepHours=0.05, AHUStatus=1):  # default time step is 1/20th of an hour, or 3 minutes
        roomTemp = self.roomTemp.getSensorReading()
        if roomTemp > self.coolingStpt:  # should run through the cooling sequence
            heatingPosition = 0.
            cfmSetpoint = self.coolingOutput.calculateCoolingAirflowStpt(roomTemp - self.coolingStpt)
            self.needCoolingStatus.setActualNeedCoolingStatus(1)
            self.needHeatingStatus.setActualNeedHeatingStatus(0 + AHUStatus)
            # print "Room Temp ({0}) is greater than the cooling setpoint ({1}). Setting the airflow setpoint to {2} cfm, and the heating valve to {3}%.".format(roomTemp, self.coolingStpt, cfmSetpoint, heatingPosition)
        elif roomTemp < self.heatingStpt:  # should run through the heating sequence
            heatingPosition = self.heatingOutput.calculateHeatingOutput(self.heatingStpt - roomTemp)
            if self.unitType in ["VAVR", "VAVRE"]:
                cfmSetpoint = self.maxHeatingCFM
                self.needHeatingStatus.setActualNeedHeatingStatus(1)
                self.needCoolingStatus.setActualNeedCoolingStatus(0)
            else:
                cfmSetpoint = 0 #self.minHeatingCFM
                self.needHeatingStatus.setActualNeedHeatingStatus(0 + AHUStatus)
                self.needCoolingStatus.setActualNeedCoolingStatus(0 + AHUStatus)
            # print "Room Temp ({0}) is less than the heating setpoint ({1}). Setting the airflow setpoint to {2} cfm, and the heating valve to {3}%.".format(roomTemp, self.coolingStpt, cfmSetpoint, heatingPosition)
        else:  # are good, so do the minimum cooling / heating
            heatingPosition = 0.
            cfmSetpoint = self.minCoolingCFM
            self.needHeatingStatus.setActualNeedHeatingStatus(0 + AHUStatus)
            self.needCoolingStatus.setActualNeedCoolingStatus(0 + AHUStatus)
            # print "Room Temp ({0}) is between cooling and heating setpoints. Setting the airflow setpoint to {2} cfm, and the heating valve to {3}%.".format(roomTemp, self.coolingStpt, cfmSetpoint, heatingPosition)

        # incorporate the fan into the heating somehow? Maybe make a heater / fan class like for damper assembly?
        self.damperAssembly.updateAirflowStpt(cfmSetpoint)
        self.damperAssembly.moveDamperToMeetStpt(timeStepHours)
        self.heaterFanPackage.commandHeaterFanPackage(heatingPosition, self.damperAssembly.getActualAirflow())
        self.roomCoolingStpt.setActualRoomTempStpt(self.coolingStpt)
        self.roomHeatingStpt.setActualRoomTempStpt(self.heatingStpt)
        self.setNeedHeat()
        self.setNeedCool()

    def runAirflowSimulation(self):
        # returns a tuple of the airflow and the temperature from the terminal unit

        # update the airflow based off of the AHU static pressure, which calculates from the airflow and damper position
        boxPressureDrop = self.AHU.Duct_Static_Pressure_BeforePIU - (
                    self.AHU.Duct_Static_Pressure_AfterPIU * ((self.damperPos.getActualPosition()) / 100.))
        self.damperAssembly.setStaticPressure(boxPressureDrop)
        actualAirflow = max(0.1,
                            self.damperAssembly.getActualAirflow())  # provide a minimum, so that does not divide by 0
        actualTemp = self.AHU.AHU_DAT.getActualTemp()
        plenumTemp = self.currentFloor.plenumTemperature

        airflowAfterHeater, airTempAfterHeater= self.heaterFanPackage.calculateTUDATandAirflow(actualTemp, actualAirflow, plenumTemp, self.roomServed.designHeatingLoads)

        self.dischargeTemp.setActualTemp(airTempAfterHeater)
        self.hotWaterReturnTemperatureandFlow()

        return airflowAfterHeater, airTempAfterHeater

    def hotWaterReturnTemperatureandFlow(self):
        return self.heaterFanPackage.calculateHeaterValveOutput()

    def calculateEnergyUsage(self):
        return self.heaterFanPackage.calculateElectricalPower(self.AHU.designDuctPressureBeforePIU,
                                                              self.AHU.designDuctPressureAfterPIU,
                                                              self.maxHeatingCFM * 0.75)

class BuildingArea:
    def __init__(self, parentBuilding, areaName="New Area", targetFloorArea=30000, isRoof=False, randomizeSensors=True, setUpWorking=None):
        self.areaName = areaName
        self.randomizeSensors = randomizeSensors
        self.setUpWorking = setUpWorking
        self.terminalUnitList = []
        self.roomList = []
        self.parentBuilding = parentBuilding
        self.associateAHU()
        self.averageAreaTemp = 70.
        self.totalRoomVolume = 0.
        self.isRoof = (isRoof == True)

        # set up dimensions
        self.plenumHeight = 3
        self.floorArea = 0
        self.createRooms(targetFloorArea, isRoof)
        self.recalculateAreaCharacteristics()

        self.resetSimulation()

    def addTerminalUnit(self, boxRepresentation):
        self.terminalUnitList.append(boxRepresentation)

    def addRoom(self, roomRepresentation):
        self.roomList.append(roomRepresentation)
        self.floorArea += roomRepresentation.floorSize
        self.totalRoomVolume += roomRepresentation.roomVolume
        if roomRepresentation.terminalUnit is not None:
            self.addTerminalUnit(roomRepresentation.terminalUnit)

    def recalculateAreaCharacteristics(self, plenumHeight=None):
        if plenumHeight is not None:
            self.plenumHeight = plenumHeight
        self.plenumVolume = self.plenumHeight * self.floorArea
        self.plenumHeatInertia = self.plenumVolume * 0.075 * 0.24  # convert the plenum volume to pounds (0.075 lbs / ft^3), and then to Btu's to increase 1 degree (0.24 degF / (BTU * lb))

        self.roomHeatInertia = self.totalRoomVolume * 0.075 * 0.24  # convert the floor volume to pounds (0.075 lbs / ft^3), and then to Btu's to increase 1 degree (0.24 degF / (BTU * lb))

    def createRooms(self, targetFloorArea=0., isRoof=None):
        if isRoof is None:
            roofUsed = self.isRoof
        else:
            roofUsed = (isRoof == True)
        nextRoomNumber = len(self.roomList) + 1
        while self.floorArea < targetFloorArea:
            newRoom = RoomRepresentation(self.areaName + "-R" + str(nextRoomNumber),
                                         self.parentBuilding,
                                         self,
                                         self.AHU,
                                         roofUsed,
                                         randomizeSensors=self.randomizeSensors,
                                         setUpWorking=self.setUpWorking)
            self.addRoom(newRoom)
            nextRoomNumber += 1
        self.recalculateAreaCharacteristics()
        self.AHU.sizeAHU()

    def associateAHU(self, AHU_Representation=None):
        if AHU_Representation is None:
            self.AHU = AirHandler("AHU-" + self.areaName, self.parentBuilding, self, randomizeSensors=self.randomizeSensors, setUpWorking=self.setUpWorking)
        else:
            self.AHU = AHU_Representation

    def resetSimulation(self):
        self.plenumTemperature = 75.
        self.averageAreaTemp = 70.
        self.AHU.resetSimulation()
        for currRoom in self.roomList:
            currRoom.resetSimulation()

    def resetSimulationTime(self, hourOfDay=0.):
        # reset all references to time within the simulation to a new time
        # hourOfDay is a float between 0 and 24, representing the time since midnight (in hours)
        for currRoom in self.roomList:
            currRoom.resetSimulationTime(hourOfDay)

    def logConditions(self, timestamp, logAsJSON=True):
        isJSONFormat = (logAsJSON == True)
        logString = ""
        # log all of the terminal unit information
        logString = addLogToRunningList(logString, logDeviceList(timestamp, self.terminalUnitList, isJSONFormat,
                                                                 logIssueStatus=True), isJSONFormat)

        # log all of the room information
        logString = addLogToRunningList(logString, logDeviceList(timestamp, self.roomList, isJSONFormat,
                                                                 logIssueStatus=True), isJSONFormat)

        # log the air handling unit information
        logString = addLogToRunningList(logString, logDevice(timestamp, self.AHU, isJSONFormat, logIssueStatus=True),
                                        isJSONFormat)

        return logString

    def createIssueSummary(self):
        logString = ""
        # log all of the terminal unit issue summaries
        for currUnit in self.terminalUnitList:
            logString = addLogToRunningList(logString, currUnit.createIssueSummary(), logAsJSON=True)

        # add in the AHU as well
        logString = addLogToRunningList(logString, self.AHU.createIssueSummary(), logAsJSON=True)

        return logString

    def simulateArea(self, timeStepHours=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        # run the simulation for the floor, and also calculate the plenum temp.
        plenumChangeInEnergy = 0.
        averageRoomTempChange = 0.
        for currRoom in self.roomList:
            airflowFromRoom, roomTemp = currRoom.simulateRoom(timeStepHours)
            plenumChangeInEnergy += (airflowFromRoom * 1.08 * (
                        roomTemp - self.plenumTemperature) + currRoom.getHeatFlowIntoPlenum()) * timeStepHours  # convert from BTUh to BTUs
            averageRoomTempChange += (roomTemp - self.averageAreaTemp) * currRoom.roomVolume / self.totalRoomVolume
        self.AHU.simulateAHU(timeStepHours)
        self.AHU.calculateFanTemp(timeStepHours)
        self.AHU.calculateFloorEnergyUsage(timeStepHours)
        self.plenumTemperature += plenumChangeInEnergy / self.plenumHeatInertia
        self.averageAreaTemp += averageRoomTempChange
        self.AHU.AHU_RAT.setActualTemp(self.plenumTemperature)  # update the return temperature with the correct heat
        self.AHU.getOutsideTempAndHumidity(self.parentBuilding.outsideConditions)

class BuildingModel:
    def __init__(self, buildingNumber=0, targetSquareFootage=500000, randomizeSensors=False, setUpWorking=True):
        self.buildingNumber = buildingNumber
        self.randomizeSensors = randomizeSensors
        self.setUpWorking = setUpWorking
        self.buildingName = "B" + str(buildingNumber)
        self.buildingSquareFootage = 0
        self.currentTime = datetime.now()  # initialize the current time object
        self.constructionMaterialParameters = ConstructionMaterialValues()  # instance of ConstructionMaterialValues which contains construction information
        self.outsideConditions = OutsideConditionsHandler()
        self.floorList = []
        self.defaultTUHeatingType = (selectRandomElement(("PIU", "VAVR", "VAVRE", "SPIU"), (0.4, 0.2, 0.2, 0.2)), selectRandomElement(("VAV", "PIU", "VAVRE", "VAVR", "SPIU"), (0.4, 0.15, 0.15, 0.15, 0.15)))

        self.heatingOutsideParameters = OutsideConditionsHandler()  # instance of OutsideConditionsHandler which contains outside conditions
        self.heatingOutsideParameters.setDesignHeatingConditions()
        self.coolingOutsideParameters = OutsideConditionsHandler()  # instance of OutsideConditionsHandler which contains outside conditions
        self.coolingOutsideParameters.setDesignCoolingConditions()

        self.occupancySchedule = OccupancySchedule()  # instance of OccupancySchedule which contains information about the occupancy of the building at various times of the day
        self.roomTempStpts = RoomTemperatureSetpoints()
        self.setUpdateTimes()
        self.chilledWaterSystem = ChilledWaterSystem(parentBuilding=self)
        self.BoilerSystem = BoilerSystem(parentBuilding=self)

        self.createRandomBuilding(targetSquareFootage)
        self.deleteLogs()  # clear out all existing logs for this building
        self.logFileName = ""
        self.logAsJSON = True


    def addFloor(self, floorName, floorplateArea, isRoof=False):
        newFloor = BuildingArea(self, floorName, floorplateArea, isRoof, randomizeSensors=self.randomizeSensors, setUpWorking=self.setUpWorking)
        self.chilledWaterSystem.addAirHandler(newFloor.AHU)
        self.floorList.append(newFloor)
        self.buildingSquareFootage += newFloor.floorArea

    def createRandomBuilding(self, targetBuildingSquareFootage=500000):
        self.buildingSquareFootage = 0
        floorPlateArea = roundToNearest(selectRandomFromRangeList(((15000., 25000.), (25000., 35000.), (35000., 50000.)), (0.4, 0.4, 0.2)), 100.)
        numberFloors = int(max(1, round(targetBuildingSquareFootage / floorPlateArea, 0)))
        if numberFloors == 1 and targetBuildingSquareFootage < 10000.:
            floorPlateArea = targetBuildingSquareFootage  # if small footprint, use one floor with the footprint of the
        for i in range(1, numberFloors + 1):
            self.addFloor(self.buildingName + "-F" + str(i), floorPlateArea, (i == numberFloors))
        self.chilledWaterSystem.sizePlant()

    def setUpdateTimes(self, simulationStepMinutes=1., logStepMinutes=1.):
        if simulationStepMinutes <= 0:  # make sure that simulation timestep is correct
            simulationStepMinutes = 1.
        self.simulationStepHours = simulationStepMinutes / 60.
        self.logStepHours = logStepMinutes / 60.
        self.simulationTimeStep = timedelta(minutes=simulationStepMinutes)
        self.logTimeStep = timedelta(minutes=logStepMinutes)

    def resetSimulation(self):
        # grab the initial time before reset
        currentSimulationTime = self.currentTime
        # resets all of the temperatures and valve positions
        for currFloor in self.floorList:
            currFloor.resetSimulation()

        # runs the simulation for a little bit to wear in the simulation
        wearInIterations = 60  # run for 60 cycles
        for i in range(wearInIterations):
            self.outsideConditions.incrementTime(self.simulationStepHours)
            self.simulateFloors(self.simulationStepHours)
            self.chilledWaterSystem.calculate()
            self.BoilerSystem.calculate()

        # after the wear-in process, reset the time to the initial time
        self.resetSimulationTime(currentSimulationTime, False)

    def resetSimulationTime(self, newSimulationTime, resetBounds = True):
        # reset all references to time within the simulation to a new time
        self.currentTime = newSimulationTime
        self.outsideConditions.sunlightData.getSunAngles(newSimulationTime)
        self.outsideConditions.setTime(self.currentTime.year, self.currentTime.month, self.currentTime.day,
                                       self.currentTime.hour, self.currentTime.minute, resetBounds)

        # relay the time of day to all floors, rooms, etc.
        hourOfDay = newSimulationTime.hour + newSimulationTime.minute / 60. + newSimulationTime.second / 3600.
        for currFloor in self.floorList:
            currFloor.resetSimulationTime(hourOfDay)

    def createLogs(self, logAsJSON=True):
        self.logAsJSON = logAsJSON != False
        if self.logAsJSON:
            # log the current sensors and actuators in the building as a JSON file
            logFile = self.buildingName + "_Logs.json"
        else:
            # log the current sensors and actuators in the building as a csv file
            logFile = self.buildingName + "_Logs.csv"
        self.logFileName = os.path.join(simulationLogDirectory, logFile)

        fh = open(self.logFileName, "w")
        if self.logAsJSON:
            fh.write("{'datapoints': [")  # to force it to write on a new line each iteration
            self.firstLineSeparator = "\n"
            self.otherLineSeparators = ",\n"
            self.isFirstLine = True
        else:
            fh.write('"Device","Time","Point","Value"')
            self.firstLineSeparator = "\n"
            self.otherLineSeparators = "\n"
            self.isFirstLine = True
        fh.close()

    def closeLogs(self):
        if self.logAsJSON:
            fh = open(self.logFileName, "a")
            fh.write("]}")  # to force it to write on a new line each iteration
            fh.close()

    def logConditions(self):
        logString = ""
        for currFloor in self.floorList:
            logString = addLogToRunningList(logString, currFloor.logConditions(self.currentTime, self.logAsJSON), self.logAsJSON)
        logString = addLogToRunningList(logString, self.chilledWaterSystem.logConditions(self.currentTime, self.logAsJSON), self.logAsJSON)

        fh = open(self.logFileName, "a")
        if self.isFirstLine:
            logTimeSeparator = self.firstLineSeparator
            self.isFirstLine = False
        else:
            logTimeSeparator = self.otherLineSeparators
        fh.write(logTimeSeparator + logString)  # to force it to write on a new line each iteration
        fh.close()

    def createIssueSummary(self):
        logString = ""
        # log all of the floor issue summaries
        for currFloor in self.floorList:
            logString = addLogToRunningList(logString, currFloor.createIssueSummary(), logAsJSON=True)
        logString = addLogToRunningList(logString, self.chilledWaterSystem.createIssueSummary(), logAsJSON=True)

        logFile = self.buildingName + "_Issue_Summary.json"
        logFile = os.path.join(simulationLogDirectory, logFile)
        fh = open(logFile, "w")
        fh.write("{'Issue List': [\n" + logString + "\n]}")  # to force it to write on a new line each iteration
        fh.close()

    def simulateFloors(self, timeStepHours=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        for currFloor in self.floorList:
            currFloor.simulateArea(timeStepHours)

    def deleteLogs(self):
        deleteSimulationLogs(self.buildingName)

    def runSimulation(self, startDate, runDays, logAsJSON=True, runSingleTime = True):
        # clear out the simulation, and run it for a while
        if runSingleTime:
            self.createLogs(logAsJSON)
        self.resetSimulationTime(startDate)
        self.resetSimulation() # resets the time back to the start time at the end of the wear in simulation
        self.simulationEndTime = self.currentTime + timedelta(days=runDays)
        self.nextLogTime = startDate

        while self.currentTime < self.simulationEndTime:
            self.outsideConditions.incrementTime(self.simulationStepHours)
            self.simulateFloors(self.simulationStepHours)
            self.chilledWaterSystem.calculate()
            self.BoilerSystem.calculate()

            if self.nextLogTime <= self.currentTime:
                self.logConditions()
                self.nextLogTime += self.logTimeStep
            self.currentTime += self.simulationTimeStep
        self.createIssueSummary()
        if runSingleTime:
            self.closeLogs()

    def runForRandomDaysThroughYear(self, numberDays=10, logAsJSON = True):
        self.createLogs(logAsJSON)
        year = 2018
        for currDayNum in range(numberDays):
            randomDate = datetime(year, int(round(selectRandomFromRange(1, 12), 0)), int(round(selectRandomFromRange(1, 30), 0)), 0)
            self.runSimulation(randomDate, 1, logAsJSON, runSingleTime=False)
            print("Done with day {0} for building {1}".format(currDayNum+1, self.buildingName))
        self.closeLogs()

class BuildingList:
    def __init__(self, numberBuildings = 1, numberDays=10, targetSquareFootage=50000):
        self.numberBuildings = numberBuildings
        self.numberDays = numberDays
        self.targetSquareFootage = targetSquareFootage
        self.numberBuildings = 0
        self.buildings = []
        self.addRandomBuilding(targetSquareFootage, numberBuildings)
        self.deleteLogs()

    def addBuilding(self, newBuilding):
        self.buildings.append(newBuilding)
        self.numberBuildings += 1

    def addRandomBuilding(self, targetSquareFootage, numberBuildings=1):
        for i in range(numberBuildings):
            self.addBuilding(BuildingModel(self.numberBuildings, targetSquareFootage))

    def runRandomNumberBuildingsForRandomDays(self, numberDays = 10, logAsJSON =True):
        for building in self.buildings:
            building.runForRandomDaysThroughYear(numberDays=numberDays, logAsJSON=logAsJSON)


    def deleteLogs(self, deleteAllInFolder=True):
        if deleteAllInFolder == True:
            deleteSimulationLogs()
        else:
            for currBuilding in self.buildings:
                currBuilding.deleteLogs()



class analogInput:
    # class which is analogous to a sensor
    # there is an input "actual" value, and an output "reading" value
    def __init__(self, sensorName, logValue=True, readingUncertainty=0., readingOffset=0.):
        self.sensorName = sensorName
        self.readingUncertainty = readingUncertainty
        self.readingOffset = readingOffset
        self.logValue = logValue  # can be True, False, or None
        self.isOverwritten = False
        self.actualValue = None
        self.readingValue = None
        self.minReading = -1000000
        self.maxReading = 1000000
        self.readingStep = 0.0001

    def setMinMax(self, minReading=None, maxReading=None, readingStep=None):
        if minReading is not None:
            self.minReading = minReading
        if maxReading is not None:
            self.maxReading = maxReading
        if readingStep is not None:
            self.readingStep = readingStep

    def overwrite(self, value):
        self.readingValue = self.readingValue = self.boundAndQuantizeValue(value)
        self.isOverwritten = True

    def getName(self):
        return self.sensorName

    def setValue(self, actualValue):
        # takes in an actual value, and then logs a "reading", which may differ from the actual or be overwritten
        self.actualValue = actualValue
        self.createReading()

    def setValues(self, actualValue, readingValue):
        self.actualValue = actualValue
        self.readingValue = self.boundAndQuantizeValue(readingValue)

    def boundAndQuantizeValue(self, currValue):
        if currValue > self.maxReading:
            return self.maxReading
        elif currValue < self.minReading:
            return self.minReading
        else:
            return roundToNearest(currValue - self.minReading, self.readingStep) + self.minReading

    def createReading(self):
        if not self.isOverwritten:
            readValue = self.readingUncertainty * selectRandomFromRange(-1, 1) + self.readingOffset + self.actualValue
            self.readingValue = self.boundAndQuantizeValue(readValue)

    def getActual(self):
        return self.actualValue

    def readValue(self):
        return self.readingValue

    def logPoint(self, logAsJSON=True):
        # value should always be a number, so don't have to wrap the value in apostrophes
        if logAsJSON == True:
            # log in a json format
            if self.logValue == True:
                return "'" + self.sensorName + "': " + str(self.readValue())
            elif self.logValue == False:
                return ""
            else:
                return "'" + self.sensorName + "': null"
        else:
            # log in a csv format
            if self.logValue == True:
                return '"' + self.sensorName + '",' + str(self.readValue())
            elif self.logValue == False:
                return ""
            else:
                return '"' + self.sensorName + '",""'

class analogOutput:
    # class which is analogous to an actuator
    # there is an input "command" value, and an output "actual" value
    def __init__(self, actuatorName, logValue=True, readingUncertainty=0., readingOffset=0.):
        self.actuatorName = actuatorName
        self.readingUncertainty = readingUncertainty
        self.readingOffset = readingOffset
        self.logValue = logValue  # can be True, False, or None
        self.isOverwritten = False
        self.actualValue = None
        self.commandValue = None
        self.minReading = -1000000
        self.maxReading = 1000000
        self.readingStep = 0.0001

    def setMinMax(self, minReading=None, maxReading=None, readingStep=None):
        if minReading is not None:
            self.minReading = minReading
        if maxReading is not None:
            self.maxReading = maxReading
        if readingStep is not None:
            self.readingStep = readingStep

    def getName(self):
        return self.actuatorName

    def overwrite(self, value):
        self.actualValue = self.boundAndQuantizeValue(value)
        self.isOverwritten = True

    def setValue(self, commandValue):
        self.commandValue = self.boundAndQuantizeValue(commandValue)
        if not self.isOverwritten:
            readingError = self.readingUncertainty * selectRandomFromRange(-1, 1) + self.readingOffset
            self.actualValue = self.boundAndQuantizeValue(readingError + self.commandValue)

    def boundAndQuantizeValue(self, currValue):
        if currValue > self.maxReading:
            return self.maxReading
        elif currValue < self.minReading:
            return self.minReading
        else:
            return roundToNearest(currValue - self.minReading, self.readingStep) + self.minReading

    def readValue(self):
        return self.commandValue

    def getActual(self):
        return self.actualValue

    def logPoint(self, logAsJSON=True):
        # value should always be a number, so don't have to wrap the value in apostrophes
        if logAsJSON == True:
            # log in a json format
            if self.logValue == True:
                return "'" + self.actuatorName + "': " + str(self.readValue())
            elif self.logValue == False:
                return ""
            else:
                return "'" + self.actuatorName + "': null"
        else:
            # log in a csv format
            if self.logValue == True:
                return '"' + self.actuatorName + '",' + str(self.readValue())
            elif self.logValue == False:
                return ""
            else:
                return '"' + self.actuatorName + '",""'

class temperatureSensor:
    def __init__(self, parentDevice=None, sensorName="Temperature", logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.setLogValue(logValue)
        self.probabilityOfSensor = 0.7

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.25,
                                            readingOffset=selectRandomFromRange(-1.5, 1.5))
        self.analogInputPoint.setMinMax(readingStep=0.25)
        self.setActualTemp(70.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        tempStatic = self.sensorName + " static"
        tempNotCalibrated = self.sensorName + " out of calibration"
        chosenIssue = selectRandomElement([tempStatic, tempNotCalibrated], [0.5, 0.5])
        self.addIssue(chosenIssue)
        if chosenIssue == tempStatic:
            self.analogInputPoint.overwrite(
                selectRandomFromRangeList(((-200., 55.), (55., 105.), (105., 255.)), (0.05, 0.9, 0.05)))
        if chosenIssue == tempNotCalibrated:
            self.analogInputPoint.readingOffset = selectRandomFromRangeList(((-20, -1.5), (1.5, 20)), (0.5, 0.5))

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getActualTemp(self):
        return self.analogInputPoint.getActual()

    def setActualTemp(self, actualTemp):
        return self.analogInputPoint.setValue(actualTemp)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class heaterOutput:
    def __init__(self, heaterDesiredBTU, parentDevice=None, outputName="Heater", isElectric=None, hasHeater=None, logValue=True,
                 setUpWorking=None):
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.logValue = logValue
        self.analogOutputPoint = None
        self.createOutput()
        self.issueList = []
        self.isElectric = isElectric
        self.designHeaterCapacity = heaterDesiredBTU
        self.setElectricStatus(self.designHeaterCapacity, isElectric)
        self.heaterCapacity = heaterDesiredBTU
        self.HWstartTemp = 180.
        self.RATdeltaT = 3.
        self.HWstartDeltaT = 3.
        self.designDTacrosscoil = 20.
        self.heaterCommandStepSize = 1 # default to having a 1% granularity
        self.probabilityOfSensor = 0.7
        if setUpWorking is None:
            self.setHeaterParameters(self.designHeaterCapacity, hasHeater=hasHeater)
        elif setUpWorking == False:
            self.setHeaterParameters(self.designHeaterCapacity, hasHeater=hasHeater, isCorrectlySized=False, isBroken=True)
        else:
            self.setHeaterParameters(self.designHeaterCapacity, hasHeater=hasHeater, isCorrectlySized=True, isBroken=False)
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogOutputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=2.5,
                                              readingOffset=selectRandomFromRange(-5., 5.))
        self.analogOutputPoint.setMinMax(readingStep=1.)

    def setElectricStatus(self, designHeatingBTU, isElectric=None):
        self.designHeaterCapacity = designHeatingBTU
        if self.isElectric:
            numberStages = 2.
            self.heaterCommandStepSize = 100. / numberStages
            self.analogOutputPoint.setMinMax(0., 100., self.heaterCommandStepSize)
            self.analogOutputPoint.readingUncertainty = 2.5  # relatively close to the actual value
            self.analogOutputPoint.readingOffset = 0.  # is either on or off; no in between
        else:
            self.analogOutputPoint.setMinMax(0., 100., 1.)
            self.analogOutputPoint.readingUncertainty = 5.  # can vary relatively far from the actual value
            self.analogOutputPoint.readingOffset = selectRandomFromRange(-10, 10)  # can be fairly off

    def setHeaterParameters(self, designHeatingBTU, hasHeater=None, hasHWV=None, isCorrectlySized=None, isBroken=None):
        if isCorrectlySized is None:  # decide if the heater is correctly sized
            isCorrectlySized = selectRandomElement([True, False], [0.9, 0.1])
        else:
            isCorrectlySized = (isCorrectlySized != False)

        if hasHeater:
            self.heaterCapacity = designHeatingBTU

            if not isCorrectlySized:
                self.addIssue(self.outputName + " Incorrectly Sized")
                self.heaterCapacity *= selectRandomFromRangeList(((0., 0.), (0.5, 0.9), (1.1, 2.0)), (0.2, 0.4, 0.4))
            else:
                self.heaterCapacity *= selectRandomFromRange(0.9, 1.1)

            if isBroken is None:  # decide if the heater is correctly sized
                isBroken = selectRandomElement([False, True], [0.9, 0.1])
            else:
                isBroken = (isBroken != False)
            if isBroken:
                self.breakOutput()
        else:
            self.heaterCapacity = 0
            self.analogOutputPoint.overwrite(0.)
            self.analogOutputPoint.logValue = None

    def calculateTempFlowAndHWRT(self, heaterDesign, actualTemp, airflowBeforeHeater):
        if self.getActualPosition() > 0.:
            totalGPM = heaterDesign / (500. * self.designDTacrosscoil)
            self.HWFlow = totalGPM * (self.getActualPosition() / 100.)
            self.HWRT = max(self.HWstartTemp - self.designDTacrosscoil, actualTemp + self.RATdeltaT)
            heaterBTU = 500. * (self.HWstartTemp - self.HWRT) * self.HWFlow
            airTempAfterHeater = actualTemp + (heaterBTU / (airflowBeforeHeater * 1.08))
            if airTempAfterHeater > self.HWstartTemp - self.HWstartDeltaT:
                airTempAfterHeater = self.HWstartTemp - self.HWstartDeltaT
                heaterBTU = airflowBeforeHeater * 1.08 * (airTempAfterHeater - actualTemp)
                self.HWRT = self.HWstartTemp - (heaterBTU / (self.HWFlow * 500.))
        else:
            airTempAfterHeater = actualTemp
            self.HWRT = self.HWstartTemp
            self.HWFlow = 0.
        return self.HWRT, self.HWFlow, airTempAfterHeater

    def breakOutput(self):
        if self.isElectric:
            self.addIssue(self.outputName + " Stuck Closed")
            self.analogOutputPoint.overwrite(0.)
        else:
            stuckOpen = self.outputName + " Stuck Open"
            stuckClosed = self.outputName + " Stuck Closed"
            stuckInPosition = self.outputName + " Stuck in position"
            issueType = selectRandomElement([stuckOpen, stuckInPosition, stuckClosed], [0.4, 0.2, 0.4])
            self.addIssue(issueType)
            if issueType == stuckClosed:
                self.analogOutputPoint.overwrite(selectRandomFromRange(0., 20.))
            elif issueType == stuckInPosition:
                self.analogOutputPoint.overwrite(selectRandomFromRange(20., 80.))
            elif issueType == stuckOpen:
                self.analogOutputPoint.overwrite(selectRandomFromRange(80., 100.))

    def reset(self):
        self.setCommandValue(0.)

    def getHeaterOutput(self):
        return self.getActualPosition() * self.heaterCapacity / 100.  # heater position goes from 0 to 100, so normalize by 100

    def getSensorReading(self):
        return self.analogOutputPoint.readValue()

    def getActualPosition(self):
        return self.analogOutputPoint.getActual()

    def setCommandValue(self, heaterCommand):
        if self.isElectric:
            return self.analogOutputPoint.setValue(roundUpToNearest(heaterCommand, self.heaterCommandStepSize))
        else:
            return self.analogOutputPoint.setValue(heaterCommand)

    def logPoint(self, logAsJSON=True):
        return self.analogOutputPoint.logPoint(logAsJSON)

class airPressureTransducer:
    def __init__(self, parentDevice=None, sensorName="Air DP", flowCoefficient=0.67, designDP=0.2, logValue=False,
                 setUpWorking=None):

        self.logValue = logValue
        self.flowCoefficient = flowCoefficient
        self.designDP = designDP
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.probabilityOfSensor = 0.7

        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()

        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.001,
                                            readingOffset=selectRandomFromRange(-0.001, 0.001))
        self.analogInputPoint.setMinMax(minReading=-10., maxReading=10., readingStep=0.0001)
        self.setActualDP(0.)  # set a default pressure

    def setPressureRange(self, minReading=0., maxReading=5.):
        # allow for a pressure range, for more versatility in sensors
        self.analogInputPoint.setMinMax(minReading=-minReading, maxReading=maxReading)

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        airDPStatic = self.sensorName + " static"
        airDPNotCalibrated = self.sensorName + " out of calibration"
        chosenIssue = selectRandomElement([airDPStatic, airDPNotCalibrated], [0.5, 0.5])
        self.addIssue(chosenIssue)
        chosenValue = self.designDP * selectRandomFromRangeList(((-3., -1.), (-1., 0.), (0., 0.), (0., 1.), (1., 3.)),
                                                                (0.15, 0.3, 0.1, 0.3, 0.15))
        if chosenIssue == airDPStatic:
            self.analogInputPoint.overwrite(chosenValue)
        if chosenIssue == airDPNotCalibrated:
            self.analogInputPoint.readingOffset = chosenValue

    def calculateDesignAirVelocity(self):
        return returnSign(self.designDP) * abs(self.designDP) ** 0.5 * 4005 * self.flowCoefficient

    def calculateSensorAirVelocity(self):
        return returnSign(self.getSensorDP()) * abs(self.getSensorDP()) ** 0.5 * 4005 * self.flowCoefficient

    def calculateActualAirVelocity(self):
        return returnSign(self.getActualDP()) * abs(self.getActualDP()) ** 0.5 * 4005 * self.flowCoefficient

    def getSensorDP(self):
        return self.analogInputPoint.readValue()

    def getActualDP(self):
        return self.analogInputPoint.getActual()

    def setActualDP(self, actualDP):
        return self.analogInputPoint.setValue(actualDP)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class damperOutput:
    def __init__(self, parentDevice=None, outputName="Damper Pos", isVariableVolume=None, logValue=True,
                 setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.analogOutputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.setDamperParameters(isVariableVolume=isVariableVolume)
        elif setUpWorking == False:
            self.setDamperParameters(isVariableVolume=isVariableVolume, isBroken=True)
        else:
            self.setDamperParameters(isVariableVolume=isVariableVolume, isBroken=False)
        self.setLogValue(logValue)
        self.probabilityOfSensor = 0.9

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogOutputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=5.,
                                              readingOffset=selectRandomFromRange(-10., 10.))
        self.analogOutputPoint.setMinMax(minReading=0., maxReading=100., readingStep=1.)

    def setDamperParameters(self, isVariableVolume=None, isBroken=None):
        if isVariableVolume is None:  # if no damper is specified is specified, decide if should have heater
            isVariableVolume = selectRandomElement([True, False], [0.9, 0.1])
        else:
            isVariableVolume = (isVariableVolume == True)  # must be true in order to have a heater

        if isVariableVolume:
            if isBroken is None:  # decide if the damper is correctly working
                isBroken = selectRandomElement([False, True], [0.9, 0.1])
            else:
                isBroken = (isBroken != False)
            if isBroken:
                self.breakOutput()
        else:
            self.analogOutputPoint.overwrite(100.)
            self.analogOutputPoint.logValue = None  # create so that all boxes have a 100% open damper - should this be true or false?

    def breakOutput(self):
        stuckOpen = self.outputName + " Stuck Open"
        stuckClosed = self.outputName + " Stuck Closed"
        stuckInPosition = self.outputName + " Stuck in position"
        issueType = selectRandomElement([stuckOpen, stuckInPosition, stuckClosed], [0.4, 0.2, 0.4])
        self.addIssue(issueType)
        if issueType == stuckClosed:
            self.analogOutputPoint.overwrite(selectRandomFromRange(0., 20.))
        elif issueType == stuckInPosition:
            self.analogOutputPoint.overwrite(selectRandomFromRange(20., 80.))
        elif issueType == stuckOpen:
            self.analogOutputPoint.overwrite(selectRandomFromRange(80., 100.))

    def getSensorReading(self):
        return self.analogOutputPoint.readValue()

    def getActualPosition(self):
        return self.analogOutputPoint.getActual()

    def setCommandValue(self, damperCommand):
        return self.analogOutputPoint.setValue(damperCommand)

    def logPoint(self, logAsJSON=True):
        return self.analogOutputPoint.logPoint(logAsJSON)

class CWVOutput:
    def __init__(self, parentDevice=None, outputName="CWV Pos", logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.analogOutputPoint = None
        self.createOutput()
        self.issueList = []
        self.setCWVParameters(setUpWorking=setUpWorking)
        self.probabilityOfSensor = 0.7
        self.designDTacrosscoil = 10.
        self.currentDP = 14.
        self.GPM = 0.

        self.setLogValue(logValue)

    def setDesignDP(self, designDP):
        self.kc = designDP / (self.designGPM ** 2)
        self.kv = self.kc / 10. # originally was 100; reduced to reduce flow at 50%

    def setCurrentDP(self, currentDP):
        self.currentDP = currentDP
        self.GPM = self.calculateCurrentGPM(currentDP)
        return self.GPM

    def setFlowFactor(self, flowFactor = 1.):
        self.GPM *= flowFactor

    def getCurrentFlow(self):
        return self.GPM


    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogOutputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(
                selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def sizeAHUcoils(self, designCoolingLoad, designHeatingLoad):
        self.designGPM = designCoolingLoad / (500 * self.designDTacrosscoil)
        print(self.designGPM)

    def calculateCurrentGPM(self, currentDP):
        return (currentDP / (self.kc + (self.kv / (max(self.getActualPosition()/100., 0.001) ** 3.)))) ** 0.5

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=0.,
                                              readingOffset=selectRandomFromRange(-1., 1.))
        self.analogOutputPoint.setMinMax(minReading=0., maxReading=100., readingStep=1.)
        self.analogOutputPoint.setValue(0.)

    def setCWVParameters(self, setUpWorking=None):
        if setUpWorking is None:  # decide if the CWV is correctly working
            setUpWorking = selectRandomElement([False, True], [0.1, 0.9])
        else:
            setUpWorking = (setUpWorking != False)
        if not setUpWorking:
            self.breakOutput()

    def breakOutput(self):
        stuckOpen = self.outputName + " Stuck Open"
        stuckClosed = self.outputName + " Stuck Closed"
        stuckInPosition = self.outputName + " Stuck in position"
        issueType = selectRandomElement([stuckOpen, stuckInPosition, stuckClosed], [0.4, 0.2, 0.4])
        self.addIssue(issueType)
        if issueType == stuckClosed:
            self.analogOutputPoint.overwrite(selectRandomFromRange(0., 20.))
        elif issueType == stuckInPosition:
            self.analogOutputPoint.overwrite(selectRandomFromRange(20., 80.))
        elif issueType == stuckOpen:
            self.analogOutputPoint.overwrite(selectRandomFromRange(80., 100.))

    def getSensorReading(self):
        return self.analogOutputPoint.readValue()

    def getActualPosition(self):
        return self.analogOutputPoint.getActual()

    def setCommandValue(self, CWVCommand):
        return self.analogOutputPoint.setValue(CWVCommand)

    def logPoint(self, logAsJSON=True):
        # return self.analogOutputPoint.logPoint(logAsJSON)

        fullLog = ""
        fullLog = addLogToRunningLog(fullLog, self.analogOutputPoint.logPoint(logAsJSON), logAsJSON)
        fullLog = addLogToRunningLog(fullLog, "'CWHV_GPM': " + str(round(self.GPM,2)), logAsJSON)
        return fullLog

    def testValve(self):
        self.setUpWorking = True
        self.sizeAHUcoils(30000,1)
        self.setDesignDP(12)
        self.setCommandValue(25)
        newDP = 12
        print("At {0} psi and {1}% open, the valve flows {2} GPM".format(round(newDP, 1), round(self.getSensorReading(), 0), round(self.setCurrentDP(newDP), 2)))
        self.setCommandValue(50)
        newDP = 12
        print("At {0} psi and {1}% open, the valve flows {2} GPM".format(round(newDP, 1),
                                                                         round(self.getSensorReading(), 0),
                                                                         round(self.setCurrentDP(newDP), 2)))
        self.setCommandValue(75)
        newDP = 12
        print("At {0} psi and {1}% open, the valve flows {2} GPM".format(round(newDP, 1),
                                                                         round(self.getSensorReading(), 0),
                                                                         round(self.setCurrentDP(newDP), 2)))
        self.setCommandValue(90)
        newDP = 12
        print("At {0} psi and {1}% open, the valve flows {2} GPM".format(round(newDP, 1),
                                                                         round(self.getSensorReading(), 0),
                                                                         round(self.setCurrentDP(newDP), 2)))
        self.setCommandValue(100)
        newDP = 12
        print("At {0} psi and {1}% open, the valve flows {2} GPM".format(round(newDP, 1),
                                                                         round(self.getSensorReading(), 0),
                                                                         round(self.setCurrentDP(newDP), 2)))



class PumpOutput:
    def __init__(self, parentDevice=None, outputName="Pump Pos", logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.analogOutputPoint = None
        self.createOutput()
        self.issueList = []
        self.setPumpParameters(setUpWorking=setUpWorking)
        self.probabilityOfSensor = 0.7

        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogOutputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(
                selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=1.,
                                              readingOffset=selectRandomFromRange(-1., 1.))
        self.analogOutputPoint.setMinMax(minReading=0., maxReading=100., readingStep=1.)
        self.analogOutputPoint.setValue(0.)

    def setPumpParameters(self, setUpWorking=None):
        if setUpWorking is None:  # decide if the CWV is correctly working
            setUpWorking = selectRandomElement([False, True], [0.1, 0.9])
        else:
            setUpWorking = (setUpWorking != False)
        if not setUpWorking:
            self.breakOutput()

    def breakOutput(self):
        stuckOpen = self.outputName + " Stuck Open"
        stuckClosed = self.outputName + " Stuck Closed"
        stuckInPosition = self.outputName + " Stuck in position"
        issueType = selectRandomElement([stuckOpen, stuckInPosition, stuckClosed], [0.4, 0.2, 0.4])
        self.addIssue(issueType)
        if issueType == stuckClosed:
            self.analogOutputPoint.overwrite(selectRandomFromRange(0., 20.))
        elif issueType == stuckInPosition:
            self.analogOutputPoint.overwrite(selectRandomFromRange(20., 80.))
        elif issueType == stuckOpen:
            self.analogOutputPoint.overwrite(selectRandomFromRange(80., 100.))

    def getSensorReading(self):
        return self.analogOutputPoint.readValue()

    def getActualPosition(self):
        return self.analogOutputPoint.getActual()

    def setCommandValue(self, PumpCommand):
        return self.analogOutputPoint.setValue(PumpCommand)

    def logPoint(self, logAsJSON=True):
        return self.analogOutputPoint.logPoint(logAsJSON)

class speedControl:
    def __init__(self, parentDevice=None, outputName="Speed", isVariableSpeed=True, logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.analogOutputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.setSpeedParameters(isVariableSpeed=isVariableSpeed)
        elif setUpWorking == False:
            self.setSpeedParameters(isVariableSpeed=isVariableSpeed, isBroken=True)
        else:
            self.setSpeedParameters(isVariableSpeed=isVariableSpeed, isBroken=False)
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogOutputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=1.,
                                              readingOffset=selectRandomFromRange(-10., 10.))
        self.analogOutputPoint.setMinMax(minReading=0., maxReading=100., readingStep=1.)

    def setSpeedParameters(self, isVariableSpeed=None, isBroken=None):
        if isVariableSpeed is None:  # if no damper is specified is specified, decide if should have heater
            isVariableSpeed = selectRandomElement([True, False], [0.9, 0.1])
        else:
            isVariableSpeed = (isVariableSpeed == True)  # must be true in order to have a heater

        if isVariableSpeed:
            if isBroken is None:  # decide if the heater is correctly sized
                isBroken = selectRandomElement([False, True], [0.9, 0.1])
            else:
                isBroken = (isBroken != False)
            if isBroken:
                self.breakOutput()
        else:
            self.analogOutputPoint.overwrite(100.)
            self.analogOutputPoint.logValue = False

    def breakOutput(self):
        stuckOpen = self.outputName + " Stuck Open"
        stuckClosed = self.outputName + " Stuck Closed"
        stuckInPosition = self.outputName + " Stuck in position"
        issueType = selectRandomElement([stuckOpen, stuckInPosition, stuckClosed], [0.4, 0.2, 0.4])
        self.addIssue(issueType)
        if issueType == stuckClosed:
            self.analogOutputPoint.overwrite(selectRandomFromRange(0., 20.))
        elif issueType == stuckInPosition:
            self.analogOutputPoint.overwrite(selectRandomFromRange(20., 80.))
        elif issueType == stuckOpen:
            self.analogOutputPoint.overwrite(selectRandomFromRange(80., 100.))

    def getSensorReading(self):
        return self.analogOutputPoint.readValue()

    def getActualSpeed(self):
        return self.analogOutputPoint.getActual()

    def setCommandValue(self, speedCommand):
        return self.analogOutputPoint.setValue(speedCommand)

    def logPoint(self, logAsJSON=True):
        return self.analogOutputPoint.logPoint(logAsJSON)

class ductAreaSetting:
    def __init__(self, designVelocity, designCFM, parentDevice=None, sensorName="Duct Area", logValue=False,
                 setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.designVelocity = designVelocity
        self.designCFM = designCFM
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        self.boxSizeOptions = [4, 6, 8, 10, 12, 14, 16]  # set up the possible diameters, in inches
        self.boxSize = self.boxSizeOptions[0]  # will be set to the actual value later
        self.piVal = 3.14159265
        if setUpWorking is None:
            self.setAreaParameters(self.designVelocity, self.designCFM)
        elif setUpWorking == False:
            self.setAreaParameters(self.designVelocity, self.designCFM, isCorrectlySized=False, isBroken=True)
        else:
            self.setAreaParameters(self.designVelocity, self.designCFM, isCorrectlySized=True, isBroken=False)
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.0,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(readingStep=0.25)

    def setAreaParameters(self, designVelocity, designCFM, isCorrectlySized=None, isBroken=None):

        exactBoxSize = self.convertFromAreaToDiameter(designCFM / designVelocity)  # convert to diameter in inches
        self.boxSize = self.boxSizeOptions[-1]
        for currBoxSize in self.boxSizeOptions:
            if (exactBoxSize + 0.25) < currBoxSize:  # size a little over the actual size, to be on the safe side
                self.boxSize = currBoxSize
                break

        if isCorrectlySized is None:  # if correct size is specified, then randomly select
            isCorrectlySized = selectRandomElement([True, False], [0.9, 0.1])
        else:
            isCorrectlySized = (isCorrectlySized != False)
        if isBroken is None:  # decide if the heater is correctly sized
            isBroken = selectRandomElement([False, True], [0.9, 0.1])
        else:
            isBroken = (isBroken != False)

        if not isCorrectlySized:
            otherBoxSizes = self.boxSizeOptions[:]  # copy the options
            otherBoxSizes.remove(self.boxSize)
            self.boxSize = selectRandomElement(otherBoxSizes)
        self.analogInputPoint.setValue(self.convertFromDiameterToArea(
            self.boxSize))  # set the actual value - if broken, the "read" value will be updated
        if isBroken:
            self.breakSensor()

    def breakSensor(self):
        incorrectValue = self.sensorName + " Incorrect Value"
        areaDefault = self.sensorName + " Default Value"
        chosenIssue = selectRandomElement([incorrectValue, areaDefault], [0.6, 0.4])
        self.addIssue(chosenIssue)
        if chosenIssue == incorrectValue:
            otherBoxSizes = self.boxSizeOptions[:]  # copy the options
            otherBoxSizes.remove(self.boxSize)
            recordedBoxSize = selectRandomElement(otherBoxSizes)
            self.analogInputPoint.overwrite(self.convertFromDiameterToArea(recordedBoxSize))
        elif chosenIssue == areaDefault:
            self.analogInputPoint.overwrite(1.)

    def convertFromDiameterToArea(self, diameterInches):
        return (diameterInches / 2. / 12.) ** 2 * self.piVal

    def convertFromAreaToDiameter(self, areaSqFt):
        return (areaSqFt / self.piVal) ** 0.5 * 2 * 12  # convert to diameter in inches

    def getProgrammedArea(self):
        return self.analogInputPoint.readValue()

    def getActualArea(self):
        return self.analogInputPoint.getActual()

    def setActualArea(self, actualArea):
        return self.analogInputPoint.setValue(actualArea)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class constantVolumeFan:
    def __init__(self, designAirflow, parentDevice=None, outputName="Fan", logValue=True, hasFan=None, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.analogOutputPoint = None
        self.createOutput()
        self.issueList = []
        self.designAirflowCapacity = designAirflow
        self.airflowCapacity = designAirflow
        self.probabilityOfSensor = 0.7
        if setUpWorking is None:
            self.setFanParameters(self.designAirflowCapacity, hasFan=hasFan)
        elif setUpWorking == False:
            self.setFanParameters(self.designAirflowCapacity, hasFan=hasFan, isCorrectlySized=False, isBroken=True)
        else:
            self.setFanParameters(self.designAirflowCapacity, hasFan=hasFan, isCorrectlySized=True, isBroken=False)
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogOutputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=10.,  # fan should do design flow, +/- 10%
                                              readingOffset=0.)  # fan should only be on or off
        self.analogOutputPoint.setMinMax(minReading=0., maxReading=100., readingStep=100.)  # will only be on or off

    def setFanParameters(self, designAirflow, hasFan=None, isCorrectlySized=None, isBroken=None):
        if isCorrectlySized is None:  # decide if the heater is correctly sized
            isCorrectlySized = selectRandomElement([True, False], [0.9, 0.1])
        else:
            isCorrectlySized = (isCorrectlySized != False)

        if hasFan:
            self.airflowCapacity = designAirflow
            if not isCorrectlySized:
                self.addIssue(self.outputName + " Incorrectly Sized")
                self.airflowCapacity *= selectRandomFromRangeList(((0., 0.), (0.5, 0.9), (1.1, 2.0)), (0.2, 0.4, 0.4))
            else:
                self.airflowCapacity *= selectRandomFromRange(0.9, 1.1)

            if isBroken is None:  # decide if the heater is correctly sized
                isBroken = selectRandomElement([False, True], [0.9, 0.1])
            else:
                isBroken = (isBroken != False)
            if isBroken:
                self.breakOutput()
        else:
            self.airflowCapacity = 0
            self.analogOutputPoint.overwrite(0.)
            self.analogOutputPoint.logValue = None

    def breakOutput(self):
        fanOffline = self.outputName + " Offline"
        fanOverridden = self.outputName + " Overridden On"
        issueType = selectRandomElement([fanOffline, fanOverridden], [0.6, 0.4])
        self.addIssue(issueType)
        if issueType == fanOffline:
            self.analogOutputPoint.overwrite(0.)
        elif issueType == fanOverridden:
            self.analogOutputPoint.overwrite(100.)

    def reset(self):
        self.setCommandValue(0.)

    def getAirflow(self):
        return self.getActualPosition() / 100. * self.airflowCapacity

    def getSensorReading(self):
        return self.analogOutputPoint.readValue()

    def getActualPosition(self):
        return self.analogOutputPoint.getActual()

    def setCommandValue(self, fanCommand):
        return self.analogOutputPoint.setValue(fanCommand)

    def logPoint(self, logAsJSON=True):
        return self.analogOutputPoint.logPoint(logAsJSON)

class variableVolumeFan:
    def __init__(self, designAirflow, parentDevice=None, outputName="Fan", logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.analogOutputPoint = None
        self.createOutput()
        self.issueList = []
        self.designAirflowCapacity = designAirflow
        self.airflowCapacity = designAirflow
        if setUpWorking is None:
            self.setFanParameters(self.designAirflowCapacity)
        elif setUpWorking == False:
            self.setFanParameters(self.designAirflowCapacity, isCorrectlySized=False, isBroken=True)
        else:
            self.setFanParameters(self.designAirflowCapacity, isCorrectlySized=True, isBroken=False)

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=5.,  # fan should be at command, +/- 5%
                                              readingOffset=0.)  # fan should only be on or off
        self.analogOutputPoint.setMinMax(minReading=0., maxReading=100., readingStep=1.)

    def setFanParameters(self, designAirflow, isCorrectlySized=None, isBroken=None):

        if isCorrectlySized is None:  # decide if the heater is correctly sized
            isCorrectlySized = selectRandomElement([True, False], [0.9, 0.1])
        else:
            isCorrectlySized = (isCorrectlySized != False)

        self.airflowCapacity = designAirflow

        if not isCorrectlySized:
            self.addIssue(self.outputName + " Incorrectly Sized")
            self.airflowCapacity *= selectRandomFromRangeList(((0., 0.), (0.5, 0.9), (1.1, 2.0)), (0.2, 0.4, 0.4))
        else:
            self.airflowCapacity *= selectRandomFromRange(0.9, 1.1)

        if isBroken is None:  # decide if the heater is correctly sized
            isBroken = selectRandomElement([False, True], [0.9, 0.1])
        else:
            isBroken = (isBroken != False)
        if isBroken:
            self.breakOutput()

    def breakOutput(self):
        fanOffline = self.outputName + " Offline"
        fanOverridden = self.outputName + " Overridden On"
        fanSetAtSpeed = self.outputName + " Set at Speed"
        issueType = selectRandomElement([fanOffline, fanOverridden], [0.6, 0.4])
        self.addIssue(issueType)
        if issueType == fanOffline:
            self.analogOutputPoint.overwrite(0.)
        elif issueType == fanOverridden:
            self.analogOutputPoint.overwrite(100.)
        elif issueType == fanSetAtSpeed:
            self.analogOutputPoint.overwrite(selectRandomFromRange(10., 90.))

    def reset(self):
        self.setCommandValue(0.)

    def getAirflow(self):
        return self.getActualPosition() / 100. * self.airflowCapacity

    def getSensorReading(self):
        return self.analogOutputPoint.readValue()

    def getActualPosition(self):
        return self.analogOutputPoint.getActual()

    def setCommandValue(self, fanCommand):
        return self.analogOutputPoint.setValue(fanCommand)

    def logPoint(self, logAsJSON=True):
        return self.analogOutputPoint.logPoint(logAsJSON)

class damperAirflow:
    def __init__(self, designAirflow, designDP, parentDevice=None, outputName="Airflow", isVariableVolume=None,
                 logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.designAirflow = designAirflow
        self.designDP = designDP
        self.ductStaticPressure = self.designDP
        self.issueList = []
        self.airflowSetpoint = 0
        self.probabilityOfSensor = 0.7
        self.staticPressureSensor = airPressureTransducer(parentDevice=self, designDP=self.designDP,
                                                          setUpWorking=setUpWorking)
        self.damperPosition = damperOutput(parentDevice=self, isVariableVolume=isVariableVolume,
                                           setUpWorking=setUpWorking, logValue=logValue)
        self.damperPosition.probabilityOfSensor = 0.9
        self.ductArea = ductAreaSetting(self.staticPressureSensor.calculateDesignAirVelocity(), self.designAirflow,
                                        parentDevice=self, setUpWorking=setUpWorking)
        self.TUairflowSP = TUairflowSP(parentDevice=self, setUpWorking=setUpWorking)
        self.analogInputPoint = None
        self.createOutput()

        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))
        relevantPoints = [self.damperPosition, self.TUairflowSP]
        for point in relevantPoints:
            point.randomizeSensorStatus()

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.staticPressureSensor.clearIssues()
        self.damperPosition.clearIssues()
        self.ductArea.clearIssues()
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.outputName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., readingStep=1.)

    def reset(self):
        self.ductStaticPressure = 0.
        self.commandDamper(0.)

    def updateAirflowStpt(self, newAirflowStpt):
        self.airflowSetpoint = newAirflowStpt
        self.TUairflowSP.setTUairflowSP(self.airflowSetpoint)

    def moveDamperToMeetStpt(self, timeStepHours=0.05):  # default time step is 1/20th of an hour, or 3 minutes
        maxStrokeMovement = 300. * timeStepHours  # assume can stroke 50% in one minute, or 3000% per hour
        currentDamperPosition = max(1.,
                                    self.damperPosition.getSensorReading())  # place a lower bound of 1 in case of a zero
        currentAirflow = max(1., self.getSensorAirflow())  # place a lower bound of 1 in case of a zero
        necessaryMovement = (self.airflowSetpoint - currentAirflow) / currentAirflow * currentDamperPosition
        if necessaryMovement < -maxStrokeMovement:
            necessaryMovement = -maxStrokeMovement
        elif necessaryMovement > maxStrokeMovement:
            necessaryMovement = maxStrokeMovement
        self.commandDamper(
            currentDamperPosition + necessaryMovement)  # the analog output will automatically place command between 0 and 100

    def commandDamper(self, damperCommand):
        self.damperPosition.setCommandValue(damperCommand)
        self.setStaticPressure()  # update the static pressure within the box

    def setStaticPressure(self, ductStaticPressure=None):
        # Actual static pressure will read correlated to the damper position.
        # Damper positions is proportional to the airflow, and airflow is proportional to velocity.
        # Static pressure is proportional to the square of the velocity
        # Therefore, the actual static pressure will be proportional to the square of the normalized damper position
        if ductStaticPressure is not None:
            self.ductStaticPressure = ductStaticPressure
            self.staticPressureSensor.setActualDP(ductStaticPressure)  # set the pressure drop across the box
        self.updateAirflows()  # update the airflows going through the box

    def updateAirflows(self):
        actualAirflow = self.staticPressureSensor.calculateActualAirVelocity() * self.ductArea.getActualArea() * self.damperPosition.getActualPosition() / 100.
        readingAirflow = self.staticPressureSensor.calculateSensorAirVelocity() * self.ductArea.getProgrammedArea() * self.damperPosition.getSensorReading() / 100.
        self.analogInputPoint.setValues(actualAirflow, readingAirflow)

    def getSensorAirflow(self):
        return self.analogInputPoint.readValue()

    def getActualAirflow(self):
        return self.analogInputPoint.getActual()

    def logPoint(self, logAsJSON=True):
        fullLog = ""
        fullLog = addLogToRunningLog(fullLog, self.damperPosition.logPoint(logAsJSON), logAsJSON)
        fullLog = addLogToRunningLog(fullLog, self.analogInputPoint.logPoint(logAsJSON), logAsJSON)
        fullLog = addLogToRunningLog(fullLog, self.ductArea.logPoint(logAsJSON), logAsJSON)
        fullLog = addLogToRunningLog(fullLog, self.staticPressureSensor.logPoint(logAsJSON), logAsJSON)
        fullLog = addLogToRunningLog(fullLog, self.TUairflowSP.logPoint(logAsJSON), logAsJSON)
        return fullLog

class heaterFanPackage:
    def __init__(self, designAirflow, heaterDesiredBTU, TUHeatingType, parentDevice=None, outputName="Heater and Fan Package", logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.TUHeatingType = TUHeatingType
        self.designAirflowCapacity = designAirflow
        self.designHeaterCapacity=heaterDesiredBTU
        self.pointsList = []
        self.pointNameList = []
        self.issueList = []
        self.minimumAirflow = 50.

        self.setUpHeater(setUpWorking=setUpWorking)
        self.setUpFan(setUpWorking=setUpWorking)  # create a fan that is either on or off

        self.analogOutputPoint = None
        self.createOutput()
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue

    def randomizeSensorStatus(self):
        for point in self.pointsList:
            point.randomizeSensorStatus()

    def getName(self):
        return self.outputName

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=2.5,
                                              readingOffset=selectRandomFromRange(-5., 5.))
        self.analogOutputPoint.setMinMax(readingStep=1.)

    def setUpHeater(self, setUpWorking=None):
        if self.TUHeatingType == "VAVR":
            self.heater = heaterOutput(self.designHeaterCapacity, parentDevice=self, outputName="Heater Pos",
                                   isElectric=False, hasHeater=True, logValue=True, setUpWorking=setUpWorking)
            self.probabilityOfSensor = 0.7
            self.addPoint(self.heater)
        elif self.TUHeatingType in ["VAVRE", "SPIU", "PIU"]:
            self.heater = heaterOutput(self.designHeaterCapacity, parentDevice=self, outputName="Heater Pos",
                                   isElectric=True, hasHeater=True, logValue=True, setUpWorking=setUpWorking)
            self.probabilityOfSensor = 0.7
            self.addPoint(self.heater)
        else:
            self.heater = heaterOutput(self.designHeaterCapacity, parentDevice=self, outputName="Heater Pos",
                                       isElectric=False, hasHeater=False, logValue=None, setUpWorking=setUpWorking)
            self.probabilityOfSensor = 0.7
            self.heater.setCommandValue(0.) # since does not have a heater, initialize to 0 so that never changes
            self.addPoint(self.heater)

    def setUpFan(self, setUpWorking=None):
        if self.TUHeatingType in ["VAV", "VAVR", "VAVRE"]:
            self.fanStatus = constantVolumeFan(0., parentDevice=self, outputName="Fan",
                                           hasFan=False, logValue=None, setUpWorking=setUpWorking)
            self.probabilityOfSensor = 0.7
            self.fanStatus.setCommandValue(0.) # since does not have a fan, initialize to 0 so that never changes
            self.addPoint(self.fanStatus)
        else:
            self.fanStatus = constantVolumeFan(self.designAirflowCapacity, parentDevice=self, outputName="Fan",
                                               hasFan=True, logValue=True, setUpWorking=setUpWorking)
            self.probabilityOfSensor = 0.7
            self.addPoint(self.fanStatus)

    def commandHeaterFanPackage(self, heatingPosition, airflow, numberofCoils=2.):
        if self.TUHeatingType in ["VAVR", "VAVRE"]:
            self.heater.setCommandValue(heatingPosition)
        elif self.TUHeatingType == "PIU":
            numberofStages = numberofCoils + 1.
            self.fanStatus.setCommandValue(heatingPosition * 100.) # if is 1% or greater, then turn on the fan to 100%
            self.heater.setCommandValue(heatingPosition * (numberofStages / numberofCoils) - (100. / numberofCoils))
        elif self.TUHeatingType == "SPIU":
            self.fanStatus.setCommandValue(100.)
            self.heater.setCommandValue(heatingPosition)

    def calculatedMixedAirflowAndTemp(self, coolingTemp, coolingAirflow, plenumTemp):
        fanAirflow = self.fanStatus.getAirflow()
        if self.TUHeatingType == "SPIU":
            mixedAirflow = max(fanAirflow, coolingAirflow)
            mixedAirTemp = (coolingAirflow * coolingTemp + (mixedAirflow - coolingAirflow) * plenumTemp) / mixedAirflow
        else:
            mixedAirflow = fanAirflow + coolingAirflow
            mixedAirTemp = (coolingAirflow * coolingTemp + fanAirflow * plenumTemp) / mixedAirflow
            # if fanAirflow > 0:
                # print coolingAirflow, coolingTemp, fanAirflow, plenumTemp, mixedAirTemp
        return mixedAirflow, mixedAirTemp

    def calculateTUDATandAirflow(self, actualTemp, airflow, plenumTemp, heaterDesign):
        airflowBeforeHeater, airTempBeforeHeater = self.calculatedMixedAirflowAndTemp(actualTemp, airflow, plenumTemp)
        if self.TUHeatingType == "VAVR":
            self.HWRT, self.HWFlow, airTempAfterHeater = self.heater.calculateTempFlowAndHWRT(heaterDesign, actualTemp, airflowBeforeHeater)
        else:
            if airflowBeforeHeater < self.minimumAirflow:
                airTempAfterHeater = airTempBeforeHeater
            else:
                heatAddedToAirstream = self.heater.getHeaterOutput()
                airTempAfterHeater = heatAddedToAirstream / (1.08 * airflowBeforeHeater) + airTempBeforeHeater
                #if airTempAfterHeater - airTempBeforeHeater < 10 and heatAddedToAirstream > 0:
                    #print round(heatAddedToAirstream,2), round(airflow,2), round(actualTemp, 2), round(airflowBeforeHeater,2), round(airTempBeforeHeater,2), round(airTempAfterHeater,2)
        self.calculateHeaterValveOutput()
        return airflowBeforeHeater, airTempAfterHeater

    def calculateHeaterValveOutput(self):
        if self.TUHeatingType == "VAVR":
            return self.HWRT, self.HWFlow

    def calculateElectricalPower(self, designDuctPressureBeforePIU, designDuctPressureAfterPIU, designAirflow):
        #calculate fan electrical power
        designDischargeStaticPressure = (designDuctPressureAfterPIU + designDuctPressureBeforePIU) / 2.
        self.fanEfficiency = 0.85
        self.designSuctionStaticPressure = -2.5
        fanBrakeHP = designAirflow * (designDischargeStaticPressure - self.designSuctionStaticPressure) / (
                    6356. * self.fanEfficiency)
        flowFactor = 0.67  # corrects for the fricitional factors when moving air, so fan needs to be sized up to correctly handle
        fanSizingSafetyFactor = 1.5

        hpStep = 0.1

        self.fanHP = roundUpToNearest(fanBrakeHP * fanSizingSafetyFactor / flowFactor, hpStep)
        self.electricalPower = ((self.fanStatus.getActualPosition() / 100.) ** 3) * self.fanHP * 2544.

        self.totalElectricalPower = (self.electricalPower + self.heater.getHeaterOutput()) / 3.41    #in Watts
        return self.totalElectricalPower

    def addPoint(self, pointInfo):
        pointName = pointInfo.getName()
        if not (pointName in self.pointNameList):
            self.pointsList.append(pointInfo)
            self.pointNameList.append(pointName)

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def logPoint(self, logAsJSON=True):
        fullLog = ""
        fullLog = addLogToRunningLog(fullLog, self.heater.logPoint(logAsJSON), logAsJSON)
        fullLog = addLogToRunningLog(fullLog, self.fanStatus.logPoint(logAsJSON), logAsJSON)
        return fullLog

class airflowStation:
    def __init__(self, parentDevice=None, sensorName="Airflow", logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=50,
                                            readingOffset=selectRandomFromRange(-100., 100.))
        self.analogInputPoint.setMinMax(minReading=0, readingStep=1.)
        self.setActualAirflow(0.)  # set a default airflow

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        tempStatic = self.sensorName + " static"
        tempNotCalibrated = self.sensorName + " out of calibration"
        chosenIssue = selectRandomElement([tempStatic, tempNotCalibrated], [0.5, 0.5])
        self.addIssue(chosenIssue)
        if chosenIssue == tempStatic:
            self.analogInputPoint.overwrite(selectRandomFromRangeList(((0., 1000.), (1000., 10000.)), (0.4, 0.6)))
        if chosenIssue == tempNotCalibrated:
            self.analogInputPoint.readingOffset = selectRandomFromRangeList(((0., 1000.), (1000., 10000.)), (0.4, 0.6))

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getActualAirflow(self):
        return self.analogInputPoint.getActual()

    def setActualAirflow(self, actualAirflow):
        return self.analogInputPoint.setValue(actualAirflow)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class waterSystemDP:
    def __init__(self, parentDevice=None, sensorName="CHWS DP", logValue=False, setUpWorking=None):

        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.probabilityOfSensor = 0.9

        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()

        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.1,
                                            readingOffset=selectRandomFromRange(-0.2, 0.2))
        self.analogInputPoint.setMinMax(minReading=-20., maxReading=50., readingStep=0.1)
        self.setActualDP(0.)  # set a default pressure

    def setPressureRange(self, minReading=-5., maxReading=20.):
        # allow for a pressure range, for more versatility in sensors
        self.analogInputPoint.setMinMax(minReading=-minReading, maxReading=maxReading)

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        waterDPStatic = self.sensorName + " static"
        waterDPNotCalibrated = self.sensorName + " out of calibration"
        chosenIssue = selectRandomElement([waterDPStatic, waterDPNotCalibrated], [0.5, 0.5])
        self.addIssue(chosenIssue)
        dpScale = 10
        chosenValue = dpScale * selectRandomFromRangeList(((-3., -1.), (-1., 0.), (0., 0.), (0., 1.), (1., 3.)),
                                                                (0.15, 0.3, 0.1, 0.3, 0.15))
        if chosenIssue == waterDPStatic:
            self.analogInputPoint.overwrite(chosenValue)
        if chosenIssue == waterDPNotCalibrated:
            self.analogInputPoint.readingOffset = chosenValue

    def getSensorDP(self):
        return self.analogInputPoint.readValue()

    def getActualDP(self):
        return self.analogInputPoint.getActual()

    def setActualDP(self, actualDP):
        return self.analogInputPoint.setValue(actualDP)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class variableVolumeHeater:#RL created this
   #Don't know what to do here
    def __init__(self, designAirflow, parentDevice=None, outputName="VAVHeater", logValue=True, setUpWorking=None):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.outputName = outputName
        self.analogOutputPoint = None
        self.createOutput()
        self.issueList = []
        self.designAirflowCapacity = designAirflow
        self.airflowCapacity = designAirflow
        if setUpWorking is None:
            self.setFanParameters(self.designAirflowCapacity)
        elif setUpWorking == False:
            self.setFanParameters(self.designAirflowCapacity, isCorrectlySized=False, isBroken=True)
        else:
            self.setFanParameters(self.designAirflowCapacity, isCorrectlySized=True, isBroken=False)

    def getName(self):
        return self.outputName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogOutputPoint = analogOutput(self.outputName,
                                              self.logValue,
                                              readingUncertainty=5.,  # fan should be at command, +/- 5%
                                              readingOffset=0.)  # fan should only be on or off
        self.analogOutputPoint.setMinMax(minReading=0., maxReading=100., readingStep=1.)

    def setFanParameters(self, designAirflow, isCorrectlySized=None, isBroken=None):

        if isCorrectlySized is None:  # decide if the heater is correctly sized
            isCorrectlySized = selectRandomElement([True, False], [0.9, 0.1])
        else:
            isCorrectlySized = (isCorrectlySized != False)

        self.airflowCapacity = designAirflow

        if not isCorrectlySized:
            self.addIssue(self.outputName + " Incorrectly Sized")
            self.airflowCapacity *= selectRandomFromRangeList(((0., 0.), (0.5, 0.9), (1.1, 2.0)), (0.2, 0.4, 0.4))
        else:
            self.airflowCapacity *= selectRandomFromRange(0.9, 1.1)

        if isBroken is None:  # decide if the heater is correctly sized
            isBroken = selectRandomElement([False, True], [0.9, 0.1])
        else:
            isBroken = (isBroken != False)
        if isBroken:
            self.breakOutput()

    def breakOutput(self):
        fanOffline = self.outputName + " Offline"
        fanOverridden = self.outputName + " Overridden On"
        fanSetAtSpeed = self.outputName + " Set at Speed"
        issueType = selectRandomElement([fanOffline, fanOverridden], [0.6, 0.4])
        self.addIssue(issueType)
        if issueType == fanOffline:
            self.analogOutputPoint.overwrite(0.)
        elif issueType == fanOverridden:
            self.analogOutputPoint.overwrite(100.)
        elif issueType == fanSetAtSpeed:
            self.analogOutputPoint.overwrite(selectRandomFromRange(10., 90.))

    def reset(self):
        self.setCommandValue(0.)

    def getAirflow(self):
        return self.getActualPosition() / 100. * self.airflowCapacity

    def getSensorReading(self):
        return self.analogOutputPoint.readValue()

    def getActualPosition(self):
        return self.analogOutputPoint.getActual()

    def setCommandValue(self, fanCommand):
        return self.analogOutputPoint.setValue(fanCommand)

    def logPoint(self, logAsJSON=True):
        return self.analogOutputPoint.logPoint(logAsJSON)


class OccupancyStatus:
    def __init__(self, parentDevice=None, sensorName="Occupancy Status", logValue=True, setUpWorking=True):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.logValue = logValue
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., maxReading=1., readingStep=1.)
        self.setActualOccupancyStatus(0.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        pass

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getActualOccupancyStatus(self):
        return self.analogInputPoint.getActual()

    def setActualOccupancyStatus(self, actualOccupancy):
        return self.analogInputPoint.setValue(actualOccupancy)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class TUairflowSP:
    def __init__(self, parentDevice=None, sensorName="TU Airflow SP", logValue=True, setUpWorking=True):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., maxReading=10000., readingStep=1.)
        self.setTUairflowSP(0.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        pass

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getTUairflowSP(self):
        return self.analogInputPoint.getActual()

    def setTUairflowSP(self, TUairflowSP):
        return self.analogInputPoint.setValue(TUairflowSP)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class StaticPressureSP:
    def __init__(self, parentDevice=None, sensorName="Static Pressure SP", logValue=True, setUpWorking=True):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., maxReading=3., readingStep=0.001)
        self.setStaticPressureSP(0.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        pass

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getStaticPressureSP(self):
        return self.analogInputPoint.getActual()

    def setStaticPressureSP(self, staticPressure):
        return self.analogInputPoint.setValue(staticPressure)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class AHUDATSP:
    def __init__(self, parentDevice=None, sensorName="AHU DAT SP", logValue=True, setUpWorking=True):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., maxReading=100., readingStep=1.)
        self.setAHUDATSP(0.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        pass

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getAHUDATSP(self):
        return self.analogInputPoint.getActual()

    def setAHUDATSP(self, AHUDATSP):
        return self.analogInputPoint.setValue(AHUDATSP)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class NeedCoolingStatus:
    def __init__(self, parentDevice=None, sensorName="Need Cooling Status", logValue=True, setUpWorking=True):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., maxReading=1., readingStep=1.)
        self.setActualNeedCoolingStatus(0.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        pass

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getActualNeedCoolingStatus(self):
        return self.analogInputPoint.getActual()

    def setActualNeedCoolingStatus(self, actualNeedCooling):
        return self.analogInputPoint.setValue(actualNeedCooling)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class NeedHeatingStatus:
    def __init__(self, parentDevice=None, sensorName="Need Heating Status", logValue=True, setUpWorking=True):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., maxReading=1., readingStep=1.)
        self.setActualNeedHeatingStatus(0.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        pass

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getActualNeedHeatingStatus(self):
        return self.analogInputPoint.getActual()

    def setActualNeedHeatingStatus(self, actualNeedHeating):
        return self.analogInputPoint.setValue(actualNeedHeating)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class RoomCoolingStpt:
    def __init__(self, parentDevice=None, sensorName="Room Temp Stpt", logValue=True, setUpWorking=True):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., maxReading=100., readingStep=0.5)
        self.setActualRoomTempStpt(0.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        pass

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getActualRoomTempStpt(self):
        return self.analogInputPoint.getActual()

    def setActualRoomTempStpt(self, coolingStpt):
        return self.analogInputPoint.setValue(coolingStpt)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)

class RoomHeatingStpt:
    def __init__(self, parentDevice=None, sensorName="Room Temp Stpt", logValue=True, setUpWorking=True):
        self.logValue = logValue
        self.parentDevice = parentDevice
        self.sensorName = sensorName
        self.analogInputPoint = None
        self.createOutput()
        self.issueList = []
        if setUpWorking is None:
            self.randomlyChooseIfBroken()
        elif setUpWorking == False:
            self.breakSensor()
        self.probabilityOfSensor = 0.7
        self.setLogValue(logValue)

    def setLogValue(self, logValue=True):
        self.logValue = logValue
        self.analogInputPoint.logValue = logValue

    def randomizeSensorStatus(self):
        if self.logValue == True:
            self.setLogValue(selectRandomElement((True, None), (self.probabilityOfSensor, 1. - self.probabilityOfSensor)))

    def getName(self):
        return self.sensorName

    def addIssue(self, issueName="Problem Detected"):
        if not (issueName in self.issueList):
            self.issueList.append(issueName)
        if self.parentDevice is not None:
            self.parentDevice.addIssue(issueName)

    def clearIssues(self):
        if self.parentDevice is not None:
            for currIssue in self.issueList:
                self.parentDevice.issueList.remove(currIssue)
        self.issueList = []
        self.createOutput()  # clear out the manifestation of the issue

    def createOutput(self):
        self.analogInputPoint = analogInput(self.sensorName,
                                            self.logValue,
                                            readingUncertainty=0.,
                                            readingOffset=0.)
        self.analogInputPoint.setMinMax(minReading=0., maxReading=100., readingStep=0.5)
        self.setActualRoomTempStpt(0.)  # set a default temperature

    def randomlyChooseIfBroken(self):
        if selectRandomElement([False, True], [0.9, 0.1]):
            self.breakSensor()

    def breakSensor(self):
        pass

    def getSensorReading(self):
        return self.analogInputPoint.readValue()

    def getActualRoomTempStpt(self):
        return self.analogInputPoint.getActual()

    def setActualRoomTempStpt(self, heatingStpt):
        return self.analogInputPoint.setValue(heatingStpt)

    def logPoint(self, logAsJSON=True):
        return self.analogInputPoint.logPoint(logAsJSON)





class PIDLoop:
    def __init__(self):
        self.pWeight = 0.
        self.iWeight = 0.
        self.dWeight = 0.
        self.lastError = 0.
        self.cumulativeIError = 0.
        self.maxCumulativeError = 1000000.

    def setLoopConstants(self, pWeight, iWeight, dWeight=0.):
        self.pWeight = pWeight
        self.iWeight = iWeight
        self.dWeight = dWeight

    def setMaxCumulativeError(self, maxCumulativeError):
        self.maxCumulativeError = maxCumulativeError

    def logError(self, currentError, timeDifference=1.):
        # calculate the p (proportional) contributions
        pAmount = currentError * self.pWeight

        # calculate the i (integral) contributions
        self.cumulativeIError += (self.lastError + currentError) / 2. * timeDifference
        if self.cumulativeIError > self.maxCumulativeError:
            self.cumulativeIError = self.maxCumulativeError
        iAmount = self.cumulativeIError * self.iWeight

        # calculate the d (derivative) contributions
        dAmount = (currentError - self.lastError) / timeDifference * self.dWeight

        # update the last state
        self.lastError = currentError

        # return the combination of the contributions
        return pAmount + iAmount + dAmount

    def resetErrors(self):
        self.lastError = 0.
        self.cumulativeIError = 0.

class AHUFanOutput:
    # a class that represents a way to convert the room temperature above the setpoint to a CWV command position
    def __init__(self):
        self.pidLoop = PIDLoop()
        self.pidLoop.setLoopConstants(50., 0.1)  # i is much lower than p
        self.pidLoop.setMaxCumulativeError(10. / self.pidLoop.iWeight)

    def calculatefanspeedOutput(self, changeinfanspeed):
        # returns a command based on 0 to 100
        pidAmount = self.pidLoop.logError(changeinfanspeed)
        if pidAmount < -20.:  # bound output between 0 and 100
            return -20.
        elif pidAmount > 20:
            return 20.
        else:
            return pidAmount

class coolingCWVOutput:
    # a class that represents a way to convert the room temperature above the setpoint to a CWV command position
    def __init__(self):
        self.pidLoop = PIDLoop()
        self.pidLoop.setLoopConstants(1., 0.001)  # i is much lower than p
        self.pidLoop.setMaxCumulativeError(10. / self.pidLoop.iWeight)

    def calculateCWVCoolingOutput(self, actualDATminusSetpoint):
        pidAmount = self.pidLoop.logError(actualDATminusSetpoint)
        if pidAmount < -15.:  # bound output between 0 and 100
            return -15.
        elif pidAmount > 15.:
            return 15.
        else:
            return pidAmount

class coolingAirflowOutput:
    # a class that represents a way to convert the room temperature above the setpoint to a damper command position
    def __init__(self, minAirflow=0., maxAirflow=0., temperatureThrottleAmount=5.):
        # will go to max after temperatureThrottleAmount degrees above setpoint
        # declare default attributes for the class
        self.minAirflow = minAirflow
        self.maxAirflow = maxAirflow
        self.airflowIntercept = 0.
        self.airflowSlope = 0.
        self.throttleTemperatureDifference = temperatureThrottleAmount

        # implement the settings
        self.pidLoop = PIDLoop()
        self.setThrottleTemperatureDifference(temperatureThrottleAmount)
        self.setMinMaxAirflow(minAirflow, maxAirflow)

    def setMinMaxAirflow(self, minAirflow, maxAirflow):
        self.minAirflow = minAirflow
        self.maxAirflow = maxAirflow
        self.airflowIntercept = minAirflow
        self.airflowSlope = (maxAirflow - minAirflow) / 100.  # scale based on output of 0 to 100.

    def setThrottleTemperatureDifference(self, throttleTemperatureDifference=5.0):
        self.throttleTemperatureDifference = throttleTemperatureDifference
        self.pidLoop.setLoopConstants(10., 2.)  # i is much lower than p
        self.pidLoop.setMaxCumulativeError(100. / self.pidLoop.iWeight)

    def calculateCoolingOutput(self, errorAmount):
        # returns a command based on 0 to 100
        # if errorAmount < 0.:
        #     # means that does not require cooling (actually has been overcooled), so no cooling is needed
        #     self.pidLoop.cumulativeIError = 0 # reset the cumulative error
        #     return 0.
        # else:
        pidAmount = self.pidLoop.logError(errorAmount)
        if pidAmount < 0.:  # bound output between 0 and 100
            return 0.
        elif pidAmount > 100.:
            return 100.
        else:
            return pidAmount

    def calculateCoolingAirflowStpt(self, errorAmount):
        return self.airflowIntercept + self.airflowSlope * self.calculateCoolingOutput(errorAmount)

class heatingOutput:
    # a class that represents a way to convert the room temperature below the setpoint to a heating command
    def __init__(self):
        # will go to max after temperatureThrottleAmount degrees below setpoint
        self.pidLoop = PIDLoop()
        self.pidLoop.setLoopConstants(1., 0.1)  # i is much lower than p
        self.pidLoop.setMaxCumulativeError(100. / self.pidLoop.iWeight)

    def calculateHeatingOutput(self, errorAmount):
        # returns a command based on 0 to 100
        pidAmount = self.pidLoop.logError(errorAmount)
        if pidAmount < 0.:  # bound output between 0 and 100
            return 0.
        elif pidAmount > 100.:
            return 100.
        else:
            return pidAmount

class coolingPumpOutput:
    # a class that represents a way to convert the room temperature above the setpoint to a CWV command position
    def __init__(self):
        self.pidLoop = PIDLoop()
        self.pidLoop.setLoopConstants(1., 0.001)  # i is much lower than p
        self.pidLoop.setMaxCumulativeError(10. / self.pidLoop.iWeight)

    def calculatePumpCoolingOutput(self, pumpError):
        pidAmount = self.pidLoop.logError(pumpError)
        if pidAmount < -15.:  # bound output between 0 and 100
            return -15.
        elif pidAmount > 15.:
            return 15.
        else:
            return pidAmount


def logDevice(timestamp, deviceModel, logAsJSON=True, logIssueStatus=True):
    # function to take in a list of analog inputs and outputs, and generate a corresponding log of the AI/O
    isJSONFormat = (logAsJSON == True)
    logString = ""
    if isJSONFormat:
        for currPoint in deviceModel.pointsList:
            logString = addLogToRunningLog(logString, currPoint.logPoint(isJSONFormat), isJSONFormat)
        if len(logString) == 0:
            return ""
        # write as a JSON format
        if logIssueStatus == True:
            if len(deviceModel.issueList) > 0:
                isBrokenStr = "true"
            else:
                isBrokenStr = "false"
            return "{'Device': '" + deviceModel.getName() + "', 'Equipment_Type': '" + deviceModel.getEquipmentType() + "', 'Time': '" + str(timestamp) + "', 'Is_Broken': " + isBrokenStr + ", " + logString + "}"
        else:
            return "{'Device': '" + deviceModel.getName() + "', 'Equipment_Type': '" + deviceModel.getEquipmentType() + "', 'Time': '" + str(timestamp) + "', " + logString + "}"
    else:
        # return as a CSV format
        deviceRecordHeader = '"' + deviceModel.getName() + '","' + str(timestamp) + '",'
        for currPoint in deviceModel.pointsList:
            pointLog = currPoint.logPoint(isJSONFormat)
            if len(pointLog) > 0:
                logString = addLogToRunningList(logString, deviceRecordHeader + pointLog, isJSONFormat)
        if len(logString) == 0:
            return ""
        return logString

def logDeviceList(timestamp, deviceList, logAsJSON=True, logIssueStatus=True):
    isJSONFormat = (logAsJSON == True)
    logString = ""
    for currDevice in deviceList:
        logString = addLogToRunningList(logString, logDevice(timestamp, currDevice, isJSONFormat, logIssueStatus),
                                        isJSONFormat)
    return logString

def addLogToRunningList(runningList, newLog, logAsJSON=True):
    if len(newLog) == 0:
        return runningList
    elif len(runningList) == 0:
        return newLog
    else:
        if logAsJSON == True:
            # put it in a JSON format
            return runningList + "," + "\n" + newLog
        else:
            # put it in a CSV format
            return runningList + "\n" + newLog

def addLogToRunningLog(runningLog, newLog, logAsJSON=True):
    if len(newLog) == 0:
        return runningLog
    elif len(runningLog) == 0:
        return newLog
    else:
        if logAsJSON == True:
            # put it in a JSON format
            return runningLog + ", " + newLog
        else:
            # put it in a CSV format
            return runningLog + "\n,," + newLog


a = BuildingModel(1, 2000)
# b = BuildingList(3, 5, 5000)
a.runSimulation(datetime(2018, 5, 22, 10), 1.0, True)
#a.runForRandomDaysThroughYear(3, True)
#b.runRandomNumberBuildingsForRandomDays()