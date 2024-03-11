/*
 * Macro to find all z-stacks and process channels 1 and 2 for 3D maxima finding
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
		if(File.isDirectory(input + File.separator + list[i])) {
			processFolder(input + File.separator + list[i]);
		}
		if((endsWith(list[i], suffix)) & (indexOf(list[i], "_peaks_") == -1)) {
			processFile(input, list[i]);
//			dryRun(input, list[i]);
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
	// work on ch1
	run("Duplicate...", "title=ch1 duplicate channels=1-1");
	run("3D Maxima Finder", "minimmum=50 radiusxy=3 radiusz=3 noise=30");
	selectWindow("peaks_ch1");
	save(input + stub + "_peaks_ch1.tif");
	close;
	selectImage("ch1");
	close;
	
	// work on ch2
	selectWindow(orig);
	run("Duplicate...", "title=ch2 duplicate channels=2-2");
	run("3D Maxima Finder", "minimmum=50 radiusxy=3 radiusz=3 noise=30");
	selectWindow("peaks_ch2");
	save(input + stub + "_peaks_ch2.tif");
	close;
	selectImage("ch2");
	close;

	selectWindow(orig);
	close;
}

function dryRun(input, file) {
	print(input + file);
}
