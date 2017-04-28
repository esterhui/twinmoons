#!/usr/bin/env python2
# Reads float32 IQ samples from stdin (interleaved I,Q)
# and applies a frequency shift (velocity) and/or a acceleration
# to the data.
# If specified, applies doppler correction given in file
# Adapted from software written to track ExoMars EDM during
# Mars landing attempt
# Stephan.Esterhuizen@jpl.nasa.gov

import os
import sys
import collections
import numpy as np
import argparse
import logging;
from collections import defaultdict
import datetime
import dateutil.parser
from astropy.time import Time
from astropy.time import TimeDelta


def main(argv=None):
    if argv is None:
        argv = sys.argv

    if args.debug:
        logging.basicConfig(stream=sys.stderr,level=logging.DEBUG)
    else:
        logging.basicConfig(stream=sys.stderr,level=logging.INFO)

    if args.fs is None:
        print("Must provide sample rate with --fs")
        return -1
                

    model=None
    if args.predicts:
        if not args.utc:
            logging.error("Must provided utc time with predicts")
            sys.exit(-1)
        model=read_esa_predicts(args.predicts)

        # If model ends before samples start,
        # error on this condition
        if model['utc'][-1] < args.utc:
            logging.error("Predicts file ends before data start time")
            logging.error("Predicts time: %s to %s"\
                    %(model['utc'][0],model['utc'][-1]))
            sys.exit(-1)

        if model['utc'][0] > args.utc:
            logging.error("Predicts file starts after data start time")
            logging.error("Predicts time: %s to %s"\
                    %(model['utc'][0],model['utc'][-1]))
            sys.exit(-1)
    elif args.dopestimate:
        model=read_doppler_estimates(args.dopestimate)

    process_data(model)

def counterrotate(data,fs,fc,fc_rate,phi):
    """
    Complex mixer, input data is mixed with a phase model phi(t)

    phi(t) = exp(-j*2*pi(v+a))

    Where: 
        v = fc*t
        a = 0.5*fc_rate*t*t

    We esentially model an object with a velocity and acceleration. The
    entire model is applied to the given complex vector

    Returns iq_cr,fc_next,phi_next

    iq_cr    : Counter rotated data
    fc_next  : Next carrier rate (Hz)
    phi_next : Next phase (radians)

    """

    logging.debug("Mixing Model is ------------")
    logging.debug("    fs      : %f Hz"%(fs)) 
    logging.debug("    phi     : %f rad"%(phi)) 
    logging.debug("    fc      : %f Hz"%(fc)) 
    logging.debug("    fc_rate : %f Hz/s"%(fc_rate)) 
    # Compute data+1 in time so we can
    # get the next phase out
    t=np.arange(0,(len(data)+1)/fs,1.0/fs)
    t=t[0:len(data)+1]

    # Compute phase due to acceleration
    acc_phase=2*np.pi*(0.5*fc_rate*t**2.0)
    # Compute phase due to velocity
    vel_phase=2*np.pi*fc*t

    # Full phase
    phase=vel_phase+acc_phase+phi

    # Mix with negative to get it to be stationary
    mixing_vector=np.exp(-1j*phase);

    # Now do complex multiply
    iq_cr=data*mixing_vector[0:-1]

    # Next phase is
    phi_next=phase[-1]
    fc_next=fc+fc_rate*t[-1]

    logging.debug("next fc is %f"%(fc_next))

    return iq_cr,fc_next,phi_next

def read_doppler_estimates(filename):
    """
    Given file of format

    hh:mm dop_estimate_hz

    will generate a 

    00:01 dopestimate1
    01:20 dopestimate2
    03:21 dopestimate3

    will create a model with doppler_hz and doppler_rate_hz
    interpolated between the timetags.

    Returns doppler model dict spaced each 1 second
    assumes first model is for second 0
    """

    with open(filename) as f:
        lines = f.readlines()

    elapsed_sec=[]
    doppler_hz=[]
    for line in lines:
        fields=line.split();
        if len(fields)!=2:
            logging.error("Expected doppler file have format: hh:mm doppler_hz");

        time=fields[0]
        minutes,sec=time.split(':')

        # Get elapsed seconds and doppler
        elapsed_sec.append(float(minutes)*60.0+float(sec));
        doppler_hz.append(float(fields[1]))

    elapsed_sec=np.array(elapsed_sec)
    doppler_hz=np.array(doppler_hz)

    acceleration_hz_s = np.diff(doppler_hz)/np.diff(elapsed_sec)

    d={}
    d['doppler_hz'] = []
    d['doppler_rate_hz_s']= []

    # Build up a model for every point
    k=0
    dop_hz=0

    if elapsed_sec[0] != 0:
        logging.error("Doppler estimate MUST start at second 0");

    for sec in range(int(elapsed_sec[-1])+1):
        if sec==elapsed_sec[k]:
            if (k>0 and round(dop_hz) != round(doppler_hz[k])):
                logging.error("Something went wrong with doppler model")

            # Done, no more data
            if k >= len(acceleration_hz_s):
                break
            dop_hz=doppler_hz[k]
            acc_hz_s=acceleration_hz_s[k]
            k=k+1

        d['doppler_hz'].append(dop_hz)
        d['doppler_rate_hz_s'].append(acc_hz_s)

        # Propagate model
        dop_hz=dop_hz+acc_hz_s

    return d

def read_esa_predicts(filename,frequency_hz=401.585625e6):
    """
    Reads ESA predicts file, returns doppler and 
    doppler-rate model vs timetag

    Format looks like this:

    Spacecraft  ID:  EDM1  
    Station ONE ID:    99
    Tropospheric model: 0  
     Header lines 

     Returns dictionary of

     d['timestr']    : Time string array
     d['utc']  : astropy.Time array
     d['doppler_hz'] : doppler in Hz array at frequency_hz
     d['doppler_rate_hz_s'] : doppler rate in Hz/s array
    """

    with open(filename) as f:
        lines = f.readlines()

    data_start=False
    lineskip=0

    m={}
    k=0
    c=299792458; # Speed of light

    for line in lines:
        if not data_start and line.find('KM') > 0:
            data_start=True
            # Dict for data to go into
            n=len(lines)-lineskip-1
            logging.debug("Reserving for %d items"%(n))
            m['timestr']=np.zeros([n],dtype='|S23')
            m['doppler_hz']=np.zeros([n],dtype=float)
            m['doppler_rate_hz_s']=np.zeros([n],dtype=float)
        elif data_start:
            d=line.split()
            if len(d) == 11:
                timestr=d[0].replace('/','-')
                timestr+='T'
                timestr+=d[1]
                timestr+=d[2]
                range_m = float(d[3])*1000;
                range_rate_m_s = float(d[4])*1000;
                range_rate_rate_m_s_s = float(d[5])
            elif len(d) == 10:
                timestr=d[0].replace('/','-')
                timestr+='T'
                timestr+=d[1]
                range_m = float(d[2])*1000;
                range_rate_m_s = float(d[3])*1000;
                range_rate_rate_m_s_s = float(d[4])
            else:
                logging.error("Line not understood %s"%(d))
                sys.exit(-1)

            m['timestr'][k]             = timestr
            m['doppler_hz'][k]          = -range_rate_m_s/c*frequency_hz
            m['doppler_rate_hz_s'][k]   = -range_rate_rate_m_s_s/c*frequency_hz
            k+=1
        else:
            lineskip+=1
            continue

    # Convert time string to timestamp
    m['utc']=Time(m['timestr'])

    logging.debug('Done reading doppler file %s'%(filename))
    logging.info('Predicts filename : %s'%(filename))
    logging.info('    Start time : %s'%(m['utc'][0]))
    logging.info('    End time   : %s'%(m['utc'][-1]))
    logging.info('    Duration   : %s days'%(m['utc'][-1]-m['utc'][0]))
    logging.info('    Start model   : %0.1f Hz, %0.3f Hz/s'\
            %(m['doppler_hz'][0],m['doppler_rate_hz_s'][0]))
    logging.info('    End model     : %0.1f Hz, %0.3f Hz/s'\
            %(m['doppler_hz'][-1],m['doppler_rate_hz_s'][-1]))

    return m

def process_data(m):
    """ 
    Process data from stdin, m is a doppler model
    as obtained from read_esa_predicts
    """

    logging.debug("Sample rate is %0.3f Hz"%(args.fs))

    # data type for binary file
    if args.swapiq:
        dt=[('imag','float32'),('real','float32')]
    else:
        dt=[('real','float32'),('imag','float32')]

    # Number of samples to work with
    fs=float(args.fs)

    # Starting phase
    phi=0;
    samples_read=0
    elapsed_sec=0
    fc=args.fc;
    fc_rate=args.fc_rate;

    if m:
        if m.has_key('utc'):
            indarray=args.utc==m['utc']
            ind=np.where(indarray==True)
            if len(ind) is not 1:
                logging.error("Couldn't find UTC time %s in doppler file"%())
                sys.exit(-1)
            ind=ind[0]
        else:
            # Starts at 0 if no UTC time given
            ind=0

    past_predicts_warning_printed=False

    if args.save_doppler:
        doppler_out_fid=open(args.save_doppler,"w")
    else:
        doppler_out_fid=None
        

    while 1:
        # Compute number of samples to read
        n=int(round((elapsed_sec+1.0)*fs-samples_read))
        data=np.fromfile(sys.stdin,dtype=dt,count=n)

        current_utc = args.utc + TimeDelta(elapsed_sec,format='sec')

        logging.debug("%s (%d): Reading %d samples"%(current_utc,elapsed_sec,n))

        if len(data)==0:
            break

        # Reserve arrays
        iq_out = np.zeros ( [ 2*len(data) ] , dtype = np.float32 )
        iq = np.zeros ( [ len(data) ] , dtype = np.complex )

        # Make complex number
        iq.real=data['real']
        iq.imag=data['imag']

        # Predicts provided?
        if m:
            # If we provide a model, use the model doppler and doppler rate
            # Note, if args.fc is given, we will apply this args.fc offset to
            # all data in predicts file. We will NOT apply a args.fc_rate offset though
            # and only use the fc_rate from the predicts file
            try:
                fc=args.fc+m['doppler_hz'][ind+elapsed_sec]
                fc_rate=m['doppler_rate_hz_s'][ind+elapsed_sec]
                if m.has_key('utc'):
                    delta_time=current_utc-m['utc'][ind+elapsed_sec]
                    # Make sure our model time and current time align
                    delta_sec=(delta_time*86400).value
                    if delta_sec > 0.1 :
                        logging.error("Model time and current time mismatch by %f seconds"%(delta_sec))
                        sys.exit(-1)

            except IndexError:
                # Just keep using fc, fc_rate and keep propagating
                if not past_predicts_warning_printed:
                    past_predicts_warning_printed=True
                    if m.has_key('utc'):
                        logging.warning("Past predicts time %s, using model %f Hz %f Hz/s"\
                                %(m['utc'][-1],fc,fc_rate))
                    else:
                        logging.warning("Past predicts time %d, using model %f Hz %f Hz/s"\
                                %(elapsed_sec-1,fc,fc_rate))

        if doppler_out_fid:
            if m.has_key['utc']:
                doppler_out_fid.write("%s %f %f\n"%(current_utc,fc,fc_rate))
            else:
                doppler_out_fid.write("%d %f %f\n"%(elapsed_sec,fc,fc_rate))

            doppler_out_fid.flush()

        # Counter rotate
        iq_cr,fc,phi=counterrotate(iq,args.fs,fc,fc_rate,phi)

        # Interleave and dump to disk
        iq_out[::2]=iq_cr.real
        iq_out[1::2]=iq_cr.imag
        iq_out.tofile(sys.stdout,format='float32');

        # Increment counters
        elapsed_sec+=1
        samples_read+=n
    
    if doppler_out_fid:
        doppler_out_fid.close()
        
if __name__ == "__main__":

    # -------- Command line interface ------------
    parser = argparse.ArgumentParser(\
            description="Reads samples from stdin, copies to stdout and writes doppler corrections to csdr pipe")

    parser.add_argument('--predicts', help='Predicts, doppler correction file. ESA GM10 file format',type=str)
    parser.add_argument('--dopestimate', help='Doppler estimate file',type=str)
    parser.add_argument('--save_doppler', help='Save applied doppler to this file',type=str)
    parser.add_argument("--fs", type=float, help="Sample rate in Hz")
    parser.add_argument("--utc", type=Time, help="UTC Time of first sample, eg '2016-10-19T06:54:07.000'",default=Time.now())
    parser.add_argument("--fc", type=float, help="Mix with this frequency (Hz)",default=0)
    parser.add_argument("--fc_rate", type=float, help="Mix with this frequency rate (Hz/s)",default=0)
    parser.add_argument("--swapiq", action='store_true', help="Swap I/Q")
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    sys.exit(main())
