#!/usr/bin/python2.7

"""
script for checking to see if a file exists
also consider: os.path.isfile(file)
"""

import sys
import os

def check(file):
	"""
	if a file exists, return 'True,' else, 'False'
	"""
#	return os.path.isfile(file)
	try:
		open(file)
		return True
	except (OSError, IOError), e:
		return False

if __name__ == "__main__":
	file = sys.argv[1]
	print check(file)
