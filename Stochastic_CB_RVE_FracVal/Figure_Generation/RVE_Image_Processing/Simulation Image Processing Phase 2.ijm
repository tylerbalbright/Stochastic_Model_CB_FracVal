dir = getDirectory("Choose a directory")
results_dir = dir + "/Results/"
File.makeDirectory(results_dir)
files = getFileList(dir)
setBatchMode(true)

for (i=0; i<files.length; i++)
{
	if (endsWith(files[i], "/"))
		continue;
	else if (startsWith(files[i], "_"))
		continue;
	else
	{
		open(dir + files[i]);
		setOption("BlackBackground", true);
		makeRectangle(50, 50, 901, 901);
		run("Crop");
        run("Make Binary");
        run("Set Scale...", "distance=852 known=40 unit=micron");
        //run("Despeckle");
		run("Set Measurements...", "area centroid perimeter fit shape feret's skewness area_fraction redirect=None decimal=3");
		run("Analyze Particles...", "display");
		saveAs("Results", results_dir + files[i] + "_Results.csv");
        selectWindow("Results");
        run("Close");
        selectWindow(files[i]);
		close();
	}
}
