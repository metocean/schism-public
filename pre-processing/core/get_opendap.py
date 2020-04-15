#!/usr/bin/env python

# Name: get_opendap.py
# Purpose: retrieve ocean model data from an OPeNDAP server
# Author: Knut-Frode Dagestad (knutfd@met.no)
# Created: Dec 2014
# Copyright: (c) MET Norway 2014
# Licence:
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


import sys
import argparse
import subprocess
import datetime
import bisect

import numpy as np
from netCDF4 import Dataset, date2num, num2date

from opendap_datasets import sources

sourceNames = sources.keys()

# Check if NCO is available
try:
    subprocess.check_output('which ncks', shell=True)
except:
    sys.exit('ERROR: NCO (netCDF Operator) not available\n'
             'Please install from http://nco.sourceforge.net/')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Native arguments
    parser.add_argument('-source', dest='source',
                        default='usn-ncom', help='OPeNDAP URL to dataset '
                        'containing ocean data, or one of keyword '
                        'identifiers from file opendap_datasets.py: '
                        + str(sourceNames))

    parser.add_argument('-f', help='Output filename (netCDF)')

    # Arguments passed directly to motu-client
    parser.add_argument('-X', type=float, help=
                        'The max longitude')

    parser.add_argument('-z', type=str, default=0.0,
                        help='The min depth (float in the '
                        'interval [0 ; Inf])')

    parser.add_argument('-Z', type=str, default=0.0,
                        help='The max depth (float in the '
                        'interval [0 ; Inf])')



    args = parser.parse_args()
    main(args)
def get_vers(t):
    t0= datetime.datetime.strptime(t,'%Y-%m-%d %H:%M:%S')
    if t0>=datetime.datetime(1992,10,2) and t0<= datetime.datetime(1995,7,31):
        V='expt_19.0'
    elif t0>=datetime.datetime(1995,8,1) and t0<= datetime.datetime(2012,4,30): 
        V='expt_19.1'
    elif t0>=datetime.datetime(2012,5,1) and t0<= datetime.datetime(2013,7,31): 
        V='expt_90.9'
    elif t0>=datetime.datetime(2013,8,1) and t0<= datetime.datetime(2014,4,4):
        V='expt_91'
    elif t0>=datetime.datetime(2014,4,4) and t0<= datetime.datetime(2016,4,18):
        V='expt_91.1'
    elif t0>=datetime.datetime(2016,4,18) and t0<= datetime.datetime(2018,11,1): 
        V='expt_91.2'
    elif t0>datetime.datetime(2018,11,1) :
        V='expt_93.0'
    else:
        print("COULD NOT FIND HYCOM DATASET")
        sys.exit(-1)

    return V


def main(f,t,T,y,Y,x,X,var='',z=0,Z=0,Source='usn-ncom'):


    try:
        url = sources[Source]
    except:
        url = Source

    source = Source

    if source=='GLBu0.08':
        Vers=get_vers(t)
        url=url+Vers
    elif source=='GLBv0.08':
        Vers='expt_93.0'
        url=url+Vers
    # Parse time

    if t is not None:
        try:  # Date and time is given
            startTime = datetime.datetime.strptime(
                t, '%Y-%m-%d %H:%M:%S')
            endTime = datetime.datetime.strptime(
                T, '%Y-%m-%d %H:%M:%S')
        except:  # Only date is given
            try:
                startTime = datetime.datetime.strptime(
                    t, '%Y-%m-%d')
                endTime = datetime.datetime.strptime(
                    T, '%Y-%m-%d')
            except:
                print('Could not parse time: ' + str(t))

    # Get metadata from URL
    d = Dataset(url, 'r')

    if f is None:
        print('==================================================')
        print(d)
    print('==================================================')

    # Get dimensions and their vectors
    for dimName in d.dimensions:
        print(dimName)
        dimvals = d.variables[dimName][:]
        if dimName == 'time':
            dimvals = num2date(dimvals, units=d.variables[dimName].units)
            if t is None:  # Return first time steps, if time not given
                startTime = dimvals[-2]
                endTime = dimvals[-1]
            if startTime > dimvals[-1] or endTime < dimvals[0]:
                sys.exit('Requested time span (%s - %s) not covered '
                         'by dataset (%s - %s)' %
                         (startTime, endTime, dimvals[0], dimvals[-1]))
            timeIndexStart = bisect.bisect_left(dimvals, startTime)
            timeIndexEnd = bisect.bisect_right(dimvals, endTime)
            timeIndexStart = np.maximum(timeIndexStart, 0)
            timeIndexEnd = np.minimum(timeIndexEnd, len(dimvals)-1)
        if dimName == 'depth':
            print(str(dimvals))
            hasDepth = True
        else:
            print('\t' + str(dimvals[0]) + ' (min)')
            print('\t' + str(dimvals[-1]) + ' (max)')




        if dimName == 'lon':
            if x is None:
                x = dimvals[0]
                X = dimvals[1]
        if dimName == 'lat':
            if y is None:
                y = dimvals[0]
                Y = dimvals[1]

    # Print available parameters
    print('\nParameters (CF standard name):')
    for varName in d.variables:
        try:
            standard_name = d.variables[varName].standard_name
            print('\t' + varName + ' (' + standard_name + ')')
        except:
            print('\t' + varName)


    if var!='':
        varCommand='-v '+','.join(var)
    else:
        varCommand=''

    d.close()

    # Suggest command, if filename is not given
    if f is not None:
        # Subsetting in time and space
        subset = varCommand+' -d lon,%.2f,%.2f -d lat,%.2f,%.2f -d time,%i,%i' % \
            (x, X, y, Y, timeIndexStart, timeIndexEnd)
        if 'hasDepth' in locals():
            z = np.float(z)
            Z = np.float(Z)
            subset = subset + ' -d depth,%.2f,%.2f ' % (z, Z)

        out = ' -o out.nc '
        ncoCommand = 'ncks ' + subset + out + url + ' --overwrite'
        ncoCommand = ncoCommand + ' -o ' + f
        print(ncoCommand)
        subprocess.call(ncoCommand, shell=True)
        #subprocess.call('ncdump out.nc -v lon|tail -3', shell=True)
        #subprocess.call('ncdump out.nc -v lat|tail -3', shell=True)
    else:
        templateCommand = '%s -source %s -t \'%s\' -T \'%s\' ' \
                          ' -x %s -X %s -y %s -Y %s -f %s' % \
                          (sys.argv[0], source, startTime, endTime,
                           x, X, y, Y, 'out.nc')
        templateCommand=templateCommand+varCommand
        if 'hasDepth' in locals():
            templateCommand = templateCommand + \
                ' -z %s -Z %s ' % (z, Z)

        print('---------------------------------------------------------')
        print('Template command for data download (cut, paste and edit):')
        print(templateCommand)
        print('=========================================================')
