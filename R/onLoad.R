# 
#  Package glmmPower calculates power for the general linear 
#  multivariate model
#  Copyright (C) 2013 Sarah Kreidler.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

.onLoad <- function(libname, pkgname) {
  library(rJava)
  library(magic)
  library(car)
  library(MASS)

  # initialize the java environment and load necessary jars
  .jpackage(pkgname, lib.loc=libname)
  
  
  # initialize the java environment
  .jpackage("glmmPower")
  .jinit()
  libdir = paste(c(getwd(), "/../inst/lib/"), collapse="")
  
  .jaddClassPath(paste(c(libdir, "commons-math3-3.2.jar"), collapse=""))
  .jaddClassPath(paste(c(libdir, "json-simple-1.1.1.jar"), collapse=""))
  .jaddClassPath(paste(c(libdir, "com.kreidles.covariatepowersimulation-1.0.0.jar"), collapse=""))
  .jaddClassPath(paste(c(libdir, "edu.ucdenver.bios.javastatistics-1.2.0.jar"), collapse=""))
  .jaddClassPath(paste(c(libdir, "jsc.jar"), collapse=""))
  
}