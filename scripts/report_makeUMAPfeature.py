#!/usr/bin/env python3
"""Len Taing 2025 (TGBTG)
Script to collect the umap feature plots, copy them into a plotting dir, and generate
a report.py compatible PLOT table, i.e. tsv with paths to plotting dir
"""

import os
import sys
import shutil
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [png file 1] -f [png file 2] ...-f [png file N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of png files to put in the report plot table")
    optparser.add_option("-r", "--report_dir", help="report dir path")
    optparser.add_option("-p", "--png_dir", help="sub-dir or report dir to copy over the png files to")
    optparser.add_option("-o", "--output", help="output .tsv file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.report_dir or not options.png_dir:
        optparser.print_help()
        sys.exit(-1)

    png_dir_path = os.path.join(options.report_dir, options.png_dir)
    if not os.path.exists(png_dir_path):
        os.mkdir(png_dir_path)
    with open(options.output, "w") as out:
        #NOTE: report.py requires at least 2 cols!
        #write a header
        out.write(f"{'\t'.join(['','Feature Plot'])}\n")
        #for each of the files, write an entry in the output img:{png_dir}/{fname}
        for f in options.files:
            fname = os.path.basename(f)
            dest = os.path.join(png_dir_path, fname)
            #copy over file
            shutil.copyfile(f, dest)
            out.write(f"-\timg:{options.png_dir}/{fname}\n")

if __name__=='__main__':
    main()
