/*
 * Macro to find all z-stacks and process channel 3 to segment nucleus and remove seeds
 */

#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".tif") suffix

setBatchMode(true);
processFolder(input);
setBatchMode(false);
// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			if(indexOf(list[i], "_peaks_") == -1) {
				processFile(input, list[i]);
			}
	}
}

function processFile(input, file) {
	open(input + file);
	
	orig = getTitle();
	getDimensions(width, height, channels, slices, frames);
	if(slices < 3) {
		close();
		return;
	}
	stub = File.nameWithoutExtension();
	print(stub);
	// work on ch3
	run("Duplicate...", "title=ch3 duplicate channels=3-3");
	selectWindow(orig);
	close();
	selectWindow("ch3");
	run("Gaussian Blur...", "sigma=2 stack");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Invert", "stack"); // needed so that nucleus is EXCLUDED
	// now open the peaks
	ch1path = input + stub + "_peaks_ch1.tif";
	open(ch1path);
	ch1 = getTitle();
	imageCalculator("AND create stack", ch1, "ch3");
	res = getTitle();
	save(input + stub + "_peaks_mod_ch1.tif");
	close;
	selectWindow(ch1);
	close;
	
	ch2path = input + stub + "_peaks_ch2.tif";
	open(ch2path);
	ch2 = getTitle();
	imageCalculator("AND create stack", ch2, "ch3");
	res = getTitle();
	save(input + stub + "_peaks_mod_ch2.tif");
	close;
	selectWindow(ch2);
	close;
	selectWindow("ch3");
	close;
}
