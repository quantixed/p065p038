/*
 * Macro to open each ER distance stack and process
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output

setBatchMode(true);
processFolder(input);
setBatchMode(false);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	if(!endsWith(input, File.separator)) input = input + File.separator;
	if(!endsWith(output, File.separator)) output = output + File.separator;
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + list[i]))
			processFolder(input + list[i]);
		if(endsWith(list[i], ".tif")) {
			if(indexOf(list[i], "edt_er_") == 0) {
				processWrapper(input, output, list[i]);
				// print(input + list[i]);
			}
		}
	}
}

function processWrapper(input, output, file) {
	for (i = 0; i < 8; i++) {
		processFile(input, output, file, 10 + i * 10);
	}

}

function processFile(input, output, file, lim) {
	open(input + file);
	print(lim + " : " + input + file);
	orig = getTitle();
	setThreshold(0.000001, lim);
	run("Convert to Mask", "background=Light black");
	rename("img");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	if(max == 0) {
		close("*");
		return 0;
	}
	run("3D Manager");
	Ext.Manager3D_Segment(128, 255);
	Ext.Manager3D_AddImage;
	Ext.Manager3D_Measure;
	csvname = replace(file,"edt_er_","stats_" + lim + "_");
	csvname = replace(csvname,".tif",".csv");
	Ext.Manager3D_SaveResult("M", output + csvname); 
	Ext.Manager3D_CloseResult("M");
	Ext.Manager3D_SelectAll;
	Ext.Manager3D_Delete;
	Ext.Manager3D_Close();
	close("*");
	run("Collect Garbage");
}

