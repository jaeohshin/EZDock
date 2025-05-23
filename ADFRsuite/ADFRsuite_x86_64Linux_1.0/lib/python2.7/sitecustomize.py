adsroot = '/data/work/EZDock/ADFRsuite/ADFRsuite_x86_64Linux_1.0'
################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2016
##
################################################################################

# specify adsroot here
import sys, os
path = os.path.join(adsroot, "CCSBpckgs")
sys.path.append(path)

from os import getenv
if getenv('ADSPYTHONPATH'):
    sys.path.insert(0, getenv('ADSPYTHONPATH'))
    
from Support.path import setSysPath
setSysPath(path)
#sys.path.insert(0,os.path.abspath('.'))
