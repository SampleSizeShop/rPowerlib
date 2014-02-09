fastEmpiricalPower(design, glh, 10,10)

J("java.lang.System")$setProperty("java.home", "/Library/Java/JavaVirtualMachines/jdk1.7.0_13.jdk/Contents/Home")
J("java.lang.System")$getProperty("java.version")
J("java.lang.System")$getProperty("java.home")
library(rJava)
.jinit()

libdir = paste(c(getwd(), "/../inst/lib/"), collapse="")
classpath = paste(c(libdir, "com.kreidles.covariatepowersimulation-1.0.0.jar:",
                    libdir, "commons-math3-3.2.jar:",
                    libdir, "json-simple-1.1.1.jar"),                  
                  collapse="")

.jaddClassPath(paste(c(libdir, "commons-math3-3.2.jar"), collapse=""))
.jaddClassPath(paste(c(libdir, "json-simple-1.1.1.jar"), collapse=""))
.jaddClassPath(paste(c(libdir, "com.kreidles.covariatepowersimulation-1.0.0.jar"), collapse=""))


J("java.lang.System")$getProperty("java.version")

print(.jclassPath())

obj=.jnew("com/kreidles/EmpiricalPowerCalculator")
startTime <- proc.time()
test = .jcall(obj, "D", "calculateEmpiricalPower", 
               toJSON(design), toJSON(glh), as.integer(1000), as.integer(1000), evalString=FALSE,
              use.true.class=TRUE)
ellapsed = proc.time() - startTime


# startTime <- proc.time()
# simulateData(design, outputDir="../data/", blockSize=1000)
# ellapsed = proc.time() - startTime

