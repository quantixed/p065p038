/*
 * Macro to open all cluster stacks and get some basic stats on the clusters
 * the aim of this is to bypass the 3D image suite which for some reason did not give results
 * for some cluster stacks.
 */

#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".tif") suffix

setBatchMode(true);
resultsPath = input + File.separator + "allResults.csv";
f = File.open(resultsPath);
print(f, "path,obj,total\n");
File.close(f);
processFolder(input, resultsPath);
setBatchMode(false);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input, resultsPath) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i], resultsPath);
		if(endsWith(list[i], suffix))
			if(indexOf(list[i], "_peaks_clust_") >= 0) {
				processFile(input, list[i], resultsPath);
			}
	}
}

function processFile(input, file, resultsPath) {
	open(input + file);
	print(input + file);
	orig = getTitle();
	// find the "parent" directory, where the results are saved
	top = replace(resultsPath, "allResults.csv", "");
	// remove this from the path of the file we are using
	path = replace(input + file, top, "");
	// here's a cludge to get rid of double slashes in the remaining path
	wrong = File.separator + File.separator;
	right = File.separator;
	path = replace(path, wrong, right);
	// now get stats from the cluster image stack
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	nClust = max;
	// threshold all signals
	setThreshold(1, 65535, "raw");
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	pxTot = (voxelCount * mean) / 255;
	// read the existing contents from a file
	contents = File.openAsString(resultsPath);
	
	f = File.open(resultsPath);
		print(f, contents);
		print(f, path + "," + nClust + "," + pxTot + "\n");
	File.close(f);
	
	selectWindow(orig);
	close;
}
