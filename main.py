"""
Copyright (c) 2022 Joe Reed
"""

import argparse
import orekit
import time
from orekit.pyhelpers import setup_orekit_curdir
from downloader import download
from requests import get

from org.hipparchus.geometry.euclidean.threed import Rotation, RotationConvention, Vector3D
from org.hipparchus.util import FastMath
from org.hipparchus.stat.descriptive import DescriptiveStatistics
from org.orekit.attitudes import InertialProvider
from org.orekit.data import DataContext
from org.orekit.frames import Frame, Transform
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.time import AbsoluteDate
from org.orekit.utils import Constants, PVCoordinatesProvider, IERSConventions

# initialize the orekit java vm
vm = orekit.initVM()

def initOrekitData():
    """
    Initialize the orekit data, download if necessary.
    """
    filePath = download('https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip')
    
    setup_orekit_curdir(filename=filePath)

def loadTle(context:DataContext, catnr:int=25544):
    """
    Load a tle as an orekit TLE object

    Args:
        catnr (int, optional): The catalog number to retrieve. Defaults to 25544 (ISS).
    """
    r = get(f"https://celestrak.com/NORAD/elements/gp.php?CATNR={catnr}&FORMAT=TLE", headers={
        "accept":"*/*",
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/99.0.4844.84 Safari/537.36'
    })
    if not r.status_code == 200:
        raise RuntimeError(f"failed to load TLE for catalog number {catnr}")
    
    data = r.content.splitlines()
    
    return TLE(data[1], data[2], context.getTimeScales().getUTC())

def checkSample(date:AbsoluteDate, pvProv:PVCoordinatesProvider, context:DataContext,
                inertialFrame:Frame=None,
                numTests:int=10, testStep:float=10, verbose:bool=False):
    """
    Test the sample data

    Args:
        date (AbsoluteDate): date at which the epoch should be converted
        pvProv (PVCoordinatesProvider): satellite coordinate provider
        context (DataContext): orekit data context
        numTests (int, optional): number of time points to test in this sample. Default is 10
        testStep (float, optional): number of seconds between test points, Default is 10
    """
    if inertialFrame == None:
        inertialFrame = context.getFrames().getGCRF()
    ecef = context.getFrames().getITRF(IERSConventions.IERS_2010, False)
    
    # this is the baseline transform
    baselineInertialToFixed = inertialFrame.getTransformTo(ecef, date)
    
    # test time points
    testPoints = map(lambda i: date.shiftedBy(float(i) * testStep), range(0, numTests))
       
    results = []
    estTime = 0
    actTime = 0
    for point in testPoints:
        # evaluate the satellite position in the inertial and fixed frames
        inertialPv = pvProv.getPVCoordinates(point, inertialFrame)
        fixedPv = pvProv.getPVCoordinates(point, ecef)
        
        # estimate the transform
        fixedTx = Transform.IDENTITY
        deltaTime = point.durationFrom(date)
        eStart = time.time()
        if deltaTime > 0:
            rotation = Rotation(Vector3D.PLUS_K, deltaTime * Constants.WGS84_EARTH_ANGULAR_VELOCITY, RotationConvention.FRAME_TRANSFORM)
            fixedTx = Transform(point, rotation)
        
        tmp = baselineInertialToFixed.transformPVCoordinates(inertialPv)
        estimatedFixedPv = fixedTx.transformPVCoordinates(tmp)
        eStop = time.time()
        
        estTime += eStop - eStart
        
        # compute the actual transform
        aStart = time.time()
        actualTx = inertialFrame.getTransformTo(ecef, point)
        actualFixedPv = actualTx.transformPVCoordinates(inertialPv)
        aStop = time.time()
        
        actTime += aStop - aStart
        
        deltaPosMag = fixedPv.getPosition().distance(estimatedFixedPv.getPosition())
        deltaVelMag = fixedPv.getVelocity().distance(estimatedFixedPv.getVelocity())
        deltaVelAng = FastMath.toDegrees(Vector3D.angle(fixedPv.getVelocity(), estimatedFixedPv.getVelocity()))
        
        deltaActPosMag = fixedPv.getPosition().distance(actualFixedPv.getPosition())
        deltaActVelMag = fixedPv.getVelocity().distance(actualFixedPv.getVelocity())
        deltaActVelAng = FastMath.toDegrees(Vector3D.angle(fixedPv.getVelocity(), estimatedFixedPv.getVelocity()))
        
        if verbose:
            print(f" Time point: {date}")
            print(f"  deltaFromEpoch: {deltaTime}s")
            print(f"  posError:       {deltaPosMag}m")
            print(f"  velMagError:    {deltaVelMag}m")
            print(f"  velAngError:    {deltaVelAng}deg")
            print(f"  check posError:       {deltaActPosMag}m")
            print(f"  check velMagError:    {deltaActVelMag}m")
            print(f"  check velAngError:    {deltaActVelAng}deg")
        
        results.append({
            "deltaTime":deltaTime,
            "posError": deltaPosMag,
            "velMagError": deltaVelMag,
            "velAngError": deltaVelAng
        })
        
    return (results, estTime, actTime)

def main():
    """
    Main function
    """
    
    parser = argparse.ArgumentParser(description="Test a quicker ECI->ECF converstion mechanism")
    parser.add_argument('-s', '--sample-points',
                        type=int,
                        default=6,
                        help='number of unique time points to test (default 6)')
    parser.add_argument('--step',
                        help='the time-step between sample points, in minutes (default 15)',
                        type=float,
                        default=15)
    parser.add_argument('-v', '--verbose',
                        help='Print intermediate results in a verbose manner',
                        default=False,
                        action='store_true')
    parser.add_argument('-i', '--inertial-frame',
                        help="The inertial frame to use (either j2000 or gcrf)",
                        choices=['j2000', 'gcrf'],
                        default='gcrf')
    parser.add_argument('-t', '--tests-per-sample',
                        help="the number of tests to run per sample (default 10)",
                        type=int,
                        default=10)
    parser.add_argument('--test-step',
                        help="the number of seconds to step forward between each test (default 10.)",
                        type=float,
                        default=10.)
    
    args = parser.parse_args()
    
    # initialize orekit data
    initOrekitData()
    
    # get data context (not really needed ... but i'm OCD)
    context = DataContext.getDefault()
    
    # get the tle from the iss
    tle = loadTle(context)
    
    # initialze SGP4 to propagate the TLE around
    teme = context.getFrames().getTEME()
    propagator = TLEPropagator.selectExtrapolator(tle, InertialProvider.of(teme), 100., teme)
    
    # create sample times
    times = map(lambda i: tle.getDate().shiftedBy(i * float(args.step) * 60), range(0, args.sample_points))
    
    # pick the inertial frame
    inertialFrame = context.getFrames().getGCRF()
    if args.inertial_frame == 'j2000':
        inertialFrame = context.getFrames().getEME2000()
    
    results = {}
    estimateTime = 0
    actualTime = 0
    for date in times:
        (data, et, at) = checkSample(date, propagator, context,
                                     verbose=args.verbose,
                                     inertialFrame=inertialFrame,
                                     numTests=args.tests_per_sample,
                                     testStep=args.test_step)
        
        estimateTime += et
        actualTime += at
        
        for d in data:
            key = str(d['deltaTime'])
            if not key in results:
                results[key] = {
                    "deltaTime": d['deltaTime'],
                    "posError": DescriptiveStatistics([d['posError']]),
                    "velMagError":  DescriptiveStatistics([d['velMagError']]),
                    "velAngError":  DescriptiveStatistics([d['velAngError']])
                }
            else:
                results[key]['posError'].addValue(d['posError'])
                results[key]['velMagError'].addValue(d['velMagError'])
                results[key]['velAngError'].addValue(d['velAngError'])
    
    for r in results.values():
        print(f"Delta time {r['deltaTime']}")
        print(f"  posErr (m):   {r['posError'].getMean()} (stddev: {r['posError'].getStandardDeviation()})")
        print(f"  velErr (m/s): {r['velMagError'].getMean()} (stddev: {r['velMagError'].getStandardDeviation()})")
        print(f"  velErr (deg): {r['velAngError'].getMean()} (stddev: {r['velAngError'].getStandardDeviation()})")
        
    print()
    print("Total time:")
    print(f"  high-fidelity computation: {actualTime} seconds")
    print(f"  estimated computation:     {estimateTime} seconds")
    print(f"Percent improvement: {100. * (actualTime - estimateTime) / actualTime}%")

if __name__ in "__main__":
    main()