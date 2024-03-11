/*
 * Macro to find all z-stacks and use the modified peak files to do 3D cluster finding
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
	// we need to process to channels
	for(i = 1; i <= 2; i ++) {
		new = "ch" + i;
		selectImage(orig);
		run("Duplicate...", "title=" + new + " duplicate channels="+ i + "-" + i);
		// find seed image
		seedpath = input + stub + "_peaks_alt_"+ new + ".tif";
		open(seedpath);
		rename("seed");
//		arg = "seeds_threshold=10 local_background=0 local_diff=0 radius_0=2 radius_1=4 radius_2=6 weigth=0.50 radius_max=10 sd_value=1 local_threshold=[Gaussian fit] seg_spot=Classical watershed volume_min=9 volume_max=800 seeds=";
		arg = "seeds_threshold=10 local_background=0 local_diff=0 radius_0=2 radius_1=4 radius_2=6 weigth=0.50 radius_max=10 sd_value=1 local_threshold=[Gaussian fit] seg_spot=Classical volume_min=9 volume_max=800 seeds=";
		arg = arg + "seed" + " spots=" + new + " radius_for_seeds=2 output=Both";
		run("3D Spot Segmentation", arg);
		Ext.Manager3D_AddImage;
		Ext.Manager3D_Measure;
		Ext.Manager3D_SaveResult("M", input + stub + "_" + new + ".csv");
		Ext.Manager3D_CloseResult("M");
		Ext.Manager3D_SelectAll;
		Ext.Manager3D_Delete;
		Ext.Manager3D_Close();
		selectWindow("Index");
		save(input + stub + "_peaks_clust_"+new+".tif");
		close;
		selectImage(orig);
		close("\\Others");
	}
	selectImage(orig);
	close;
}
