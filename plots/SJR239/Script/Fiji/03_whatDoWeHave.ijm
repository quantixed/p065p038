/*
 * Macro to find all cluster result stacks and make an easy to view version
 */

#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".tif") suffix

setBatchMode(true);
print("\\Clear");
newImage("bigstack", "8-bit ramp", 600, 600, 1);
run("glasbey");
processFolder(input);
run("Delete Slice");
setBatchMode(false);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			if(indexOf(list[i], "_peaks_clust_") >= 0) {
				processFile(input, list[i]);
			}
	}
}

function processFile(input, file) {
	open(input + file);
	path = replace(input + file, "//","/");
	print(path);
	orig = getTitle();
	stub = File.nameWithoutExtension();
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setMinAndMax(0, max);
	run("Apply LUT", "stack");
	setOption("ScaleConversions", true);
	run("8-bit");
	run("glasbey");
	run("Z Project...", "projection=[Average Intensity]");
	new = getTitle();
	run("Set Label...", "label="+orig);
	run("Concatenate...", "  title=bigstack image1=bigstack image2=" + new + " image3=[-- None --]");
	selectWindow(orig);
	close;
}
