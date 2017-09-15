# S3TBX-Snow
Processors and tools for the retrieval of snow properties within the SEOM S3-SNOW project

How to build
------------

Make sure you have **[git](https://git-scm.com/)**, 
**[JDK 1.8](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)**, and 
**[Maven 3](https://maven.apache.org/)** installed. Make sure Maven find's the JDK by setting the environment variable `JAVA_HOME` to the directory where your JDK is installed. 

Clone or fork the repository at https://github.com/bcdev/s3tbx-snow. 
```
> git clone https://github.com/bcdev/s3tbx-snow.git
> cd s3tbx-snow
```

You can update your checked-out sources from the remote repository by running 
```
> git pull --rebase
```

Incremental build with Maven:
```
> mvn package
```

Clean build:
```
> mvn clean package
```  

If you encounter test failures:
```
> mvn clean package -DskipTests=true
```

The build creates a SNAP plugin module file `target/nbm/s3tbx-snow-<version>.nbm`.

How to install and run the processor as SNAP plugin 
---------------------------------------------------

Start SNAP (Desktop UI) and find the plugin manager in the main menu at 
> **Tools / Plugins**

Then 
* select tab **Downloaded**, 
* click button **Add Files** and 
* select the plugin module file `target/nbm/s3tbx-snow-<version>.nbm`. 
* Click **Install**, 
* then **Close** and 
* restart SNAP.

Once the Snow Albedo processor is installed into SNAP it can be run from the SNAP Desktop UI's main menu at
> **Optical / Thematic Land Processing / OLCI Snow Albedo**
  
Or in batch mode using SNAP's `gpt` command-line tool found in `${SNAP_HOME}/bin`. 
```
> gpt OLCI.SnowAlbedo -h
```  

Modifying, running and debugging the processor code
---------------------------------------------------

This section explains how to run and debug the Snow Albedo processor code from a Java IDE without having to install 
the plugin into SNAP.

You will need to install
* SNAP with the Sentinel-3 Toolbox (S3TBX) from http://step.esa.int/main/download/
* IntelliJ IDEA (Community Edition) IDE from https://www.jetbrains.com/idea/download/

Start IDEA and select **File / New / Project from Existing Sources**. Select the `pom.xml` (Maven project file) in the 
source directory. Leave all default settings as they are and click **Next** until IDEA asks for the JDK. Select the 
installed JDK from above and finish the dialog.

From the main menu select **Run / Edit Configurations**. In the dialog click the **+** (add) button and select 
**JAR Application**. Then the settings are as follows:

* **Name**: SNAP Desktop
* **Path to JAR:** `${SNAP_HOME}/snap/snap/core/snap-main.jar`
* **VM options:** `-Xmx4G -Dorg.netbeans.level=INFO -Dsun.java2d.noddraw=true -Dsun.awt.nopixfmt=true -Dsun.java2d.dpiaware=false` 
* **Program arguments:** `--userdir ${S3-SNOW_HOME}/target/testdir --clusters ${S3-SNOW_HOME}/target/nbm/netbeans/s3tbx --patches ${S3-SNOW_HOME}/$/target/classes`
* **Working directory:** `${SNAP_HOME}`

where 

* `${SNAP_HOME}` must be replaced by your SNAP installation directory
* `${S3-SNOW_HOME}` must be replaced by your s3tbx-snow project directory (where the `pom.xml` is located in)