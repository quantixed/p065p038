/*
 * Macro to open all inference stacks and make two channel composites
 * Inference stacks have three values 0 = background, 1 = Mito, 2 = ER
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output

setBatchMode(true);
setForegroundColor(255, 255, 255);
processFolder(input);
setBatchMode(false);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	if(!endsWith(input, File.separator)) input = input + File.separator;
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + list[i]))
			processFolder(input + list[i]);
		if(endsWith(list[i], ".tif")) {
			processFile(input, list[i]);
//			print(input + list[i]);
		}
	}
}

function processFile(input, file) {
	open(input + file);
	print(input + file);
	orig = getTitle();
	// duplicate mito and er
	run("Duplicate...", "title=mito duplicate");
	selectImage(orig);
	run("Duplicate...", "title=er duplicate");
	selectImage(orig);
	close();
	// process mito
	selectImage("mito");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setThreshold(1, 1, "raw");
	//setThreshold(1, 1);
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black");
	// process er
	selectImage("er");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setThreshold(2, 2, "raw");
	//setThreshold(2, 2);
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black");
	// merge, save and close in outputfolder
	run("Merge Channels...", "c1=mito c2=er create");
	save(output + File.separator + file);
	close;
}
