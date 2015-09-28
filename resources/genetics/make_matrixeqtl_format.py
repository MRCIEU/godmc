#!/usr/bin/python

import os
import sys
import time
from datetime import datetime, timedelta


#
def Chunk(L, n):
    for i in xrange(0, len(L), n):
        yield L[i:i+n]


def convert_seconds(num_seconds):
    """ convert seconds to days, hours, minutes, and seconds, as appropriate"""
    sec = timedelta(seconds=num_seconds)
    d = datetime(1,1,1) + sec
    return ("%dd %dh %dm %ds" % (d.day-1, d.hour, d.minute, d.second))


class ProgressBar:

    def __init__(self, colour="green"):
        """ Create a progress bar and initalise start time of task """
        self.start_time = time.time()
        self.colours = {"":0, "black":90, "red":91, "green":92, "yellow":93, "blue":94, "purple":95, "cyan":96, "white":97}
        self.start_colour = "\033[%sm" % (self.colours.get(colour))
        self.end_colour = "\033[0m"

    def print_progress(self, current, total, additional_info=""):
        """ Call inside for loop, passing current index and total length of iterable """
        if additional_info:
            additional_info = "[%s]" % additional_info
        current += 1 # be optimistic so we finish on 100 
        percent = float(current)/float(total) * 100
        remaining_time = convert_seconds((100 - percent) * (time.time() - self.start_time)/percent)
        percent_string = "%.1f" % percent
        bar = "|%s>%s|" % ("-" * int(percent/4), " " * (25 - int(percent/4))) 
        print "\r%s%s / %s - %s%% %s %s remaining: %s %s" % (self.start_colour, current, total, percent_string, bar, additional_info, remaining_time, self.end_colour),
        sys.stdout.flush()


#
if __name__ == '__main__':
    #
    tp = sys.argv[1]
    chunksize = int(sys.argv[2])
    
    with open(str(tp) + ".raw", "r") as H:
        for nid, l in enumerate(H):
            pass
        H.seek(0)
        
        #
        # HEADER
        #
        
        Head = [ x.split("_")[0] for x in H.next()[:-1].split(" ")[6:] ] # just SNPs
        
        i = 0
        for x in Chunk(Head, chunksize):
            i += 1
            OUT = open(str(tp) + "." + str(i) + ".tab", "w")                
            OUT.write("\t".join( [ "" ] + x ) + "\n")

        print "Writing data to " + str(i) + " chunks of " + str(chunksize) + " SNPs ..."
        
        #
        # DATA
        #
        count = 0
        pb = ProgressBar("cyan")
        for record in H:
            ALN  = record[:-1].split(" ")[1]
            data = record[:-1].split(" ")[6:]
            # print_progress(count, nid, "")
            pb.print_progress(count, nid, ALN)
            count = count + 1
            # print "Writing individual ... " + count
            i = 0
            for x in Chunk(data, chunksize):
                i += 1
                OUT = open(str(tp) + "." + str(i) + ".tab", "a")
                OUT.write("\t".join( [ ALN ] + x ) + "\n")
    
