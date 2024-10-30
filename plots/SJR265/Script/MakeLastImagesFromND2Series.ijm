/*
 * Macro to find all z-stacks and process select last frame and restack into a new hyperstack
 */

#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".nd2") suffix

setBatchMode(true);
print("\\Clear");
processFolder(input);
run("Concatenate...", "all_open");
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
		}
	}
}

function processFile(input, file) {
	if(!endsWith(input, File.separator)) input = input + File.separator;
	path = input + file;
	// Open file
	series = getTotalSeries(path);
	// Open file
	s = "open=["+ path +"] autoscale color_mode=Composite rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT";
	s = s + stringBuilder(series);
	run("Bio-Formats Importer", s);
	for (i = 0; i < series; i++) {
		print(i + "\t" + file);
	}
}

function getTotalSeries(path)	{
	run("Bio-Formats Macro Extensions");
	Ext.setId(path);
	// get number of series and total images in each series
	Ext.getSeriesCount(seriesCount);
	return seriesCount;
}

function stringBuilder(no) {
	str = "";
	for (i = 1; i < no + 1; i++) {
		str = str + " series_" + i;
	}
	for (i = 1; i < no + 1; i++) {
		str = str + " c_begin_" + i + "=1 c_end_" + i + "=4 c_step_" + i + "=1 t_begin_" + i + "=11 t_end_" + i + "=11 t_step_" + i + "=1";
	}
	return str;
}

	