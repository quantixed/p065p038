/*
 * Macro to open each inference stack and segment the mitochondria into separate objects
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
			if(indexOf(list[i], "stk_") == 0) {
				processFile(input, output, list[i]);
				// print(input + list[i]);
			}
		}
	}
}

function processFile(input, output, file) {
	open(input + file);
	print(input + file);
	orig = getTitle();
	// duplicate mito channel
	run("Duplicate...", "title=mito duplicate channels=1-1");
	// close original for RAM saving
	selectImage(orig);
	close();
	// Make binary
	selectImage("mito");
	run("Make Binary", "calculate black");
	// set pixel size
	run("Set Scale...", "distance=1 known=7.6 unit=nm");
	// set voxel size
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Properties...", "channels=1 slices=" + slices + " frames=1 pixel_width=7.6000 pixel_height=7.6000 voxel_depth=80");
	
	run("3D Manager");
	Ext.Manager3D_Segment(128, 255);
	Ext.Manager3D_AddImage();
	// close mito img
	close("mito");
	// get stats
	Ext.Manager3D_Measure();
	csvname = replace(file,"stk_","mito_stats_");
	csvname = replace(csvname,".tif",".csv");
	Ext.Manager3D_SaveResult("M", output + csvname); 
	Ext.Manager3D_CloseResult("M");
	// list voxels
	Ext.Manager3D_List();
	csvname = replace(file,"stk_","mito_vx_");
	csvname = replace(csvname,".tif",".csv");
	Ext.Manager3D_SaveResult("V", output + csvname); 
	Ext.Manager3D_CloseResult("V");
	Ext.Manager3D_SelectAll();
	Ext.Manager3D_Delete();
	Ext.Manager3D_Close();
	// close
	close("*");
	run("Collect Garbage");
}

