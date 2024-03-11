/*
 * The goal here is to process the overlap between ch1 and ch2
 */

#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".tif") suffix

setBatchMode(true);
resultsPath = input + File.separator + "allResults2.csv";
f = File.open(resultsPath);
print(f, "path,obj1,total1,obj2,total2,ch12,ch21\n");
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
			if(indexOf(list[i], "_peaks_clust_ch1") >= 0) {
				processFile(input, list[i], resultsPath);
			}
	}
}

function processFile(input, file, resultsPath) {
	open(input + file);
	print(input + file);
	ch1 = getTitle();
	open(input + replace(file,"ch1","ch2"));
	ch2 = getTitle();
	// find the "parent" directory, where the results are saved
	top = replace(resultsPath, "allResults.csv", "");
	// remove this from the path of the file we are using
	path = replace(input + file, top, "");
	// here's a cludge to get rid of double slashes in the remaining path
	wrong = File.separator + File.separator;
	right = File.separator;
	path = replace(path, wrong, right);
	// now get stats from the cluster image stack
	selectWindow(ch1);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	nClust1 = max;
	// threshold all signals
	setThreshold(1, 65535, "raw");
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black create");
	rename("mask1");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	pxTot1 = (voxelCount * mean) / 255;
	// now get stats from the other cluster image stack
	selectWindow(ch2);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	nClust2 = max;
	// threshold all signals
	setThreshold(1, 65535, "raw");
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black create");
	rename("mask2");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	pxTot2 = (voxelCount * mean) / 255;
	// next we need the overlap
	// ch1 that is also ch2
	selectWindow("mask2");
	run("Divide...", "value=255");
	imageCalculator("Multiply create stack", ch1, "mask2");
	rename("result1");
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj);
	ch12 = nb_obj;
	Ext.Manager3D_SelectAll();
	Ext.Manager3D_Delete();
	// ch2 that is also ch1
	selectWindow("mask1");
	run("Divide...", "value=255");
	imageCalculator("Multiply create stack", ch2, "mask1");
	rename("result2");
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj);
	ch21 = nb_obj;
	Ext.Manager3D_SelectAll();
	Ext.Manager3D_Delete();
	Ext.Manager3D_Close();
	
	// read the existing contents from a file
	contents = File.openAsString(resultsPath);
	
	f = File.open(resultsPath);
		print(f, contents);
		print(f, path + "," + nClust1 + "," + pxTot1 + "," + nClust2 + "," + pxTot2 + "," + ch12 + "," + ch21 + "\n");
	File.close(f);
	
	close("*");
}
