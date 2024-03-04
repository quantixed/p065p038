/*
 * Macro to find all inference stacks from nnU-Net and generate a mitochondria-ER distance map from them
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output

setBatchMode(true);
close("*");
print("\\Clear");
processFolder(input);
setBatchMode(false);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], ".tif"))
			if(indexOf(list[i], "stk_") != -1) {
				processFile(input, output, list[i]);
			}
	}
}

function processFile(input, output, file) {
	// open image
	open(input + File.separator + file);
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
	// Generate distance map and write out progress
	print("Start generating EDT");
	run("3D Distance Map", "map=EDT image=mito mask=Same threshold=0 inverse");
	print("Finished generating EDT");
	
	// save
	saveAs("Tiff", output + File.separator + replace(file, "stk_", "edt_"));
	rename("edt");
	
	// close the mito image
	selectImage("mito");
	close();
	
	// open original again
	open(input + File.separator + file);
	orig = getTitle();
	selectImage(orig);
	// duplicate ER channel
	run("Duplicate...", "title=er duplicate channels=2-2");
	// close original
	selectImage(orig);
	close();
	
	// work on ER image generate 1 or 0 image
	selectImage("er");
	run("Divide...", "value=255 stack");
	// generate a version where EDT values are transferred to the ER mask only
	imageCalculator("Multiply create 32-bit stack", "edt","er");
	// save
	saveAs("Tiff", output + File.separator + replace(file, "stk_", "edt_er_"));
	// close all
	close("*");
}