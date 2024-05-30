/*
 * This script will extract line profiles
 * Open image, store line profiles in the ROI manager
 */

lw = 3;
Dialog.create("Settings");
Dialog.addMessage("For this experiment:");
Dialog.addNumber("Line width (px)", lw);
Dialog.show();

lw = Dialog.getNumber();

retrieveCoatingData(lw);

function retrieveCoatingData(lw) {
	// location of outputs is the same as the image
	dir = getInfo("image.directory");
	fName = getInfo("image.filename");
	fPath = dir + File.separator + fName;
	resultBaseName = File.getNameWithoutExtension(fPath);
	fPath = dir + File.separator + resultBaseName; // now set this as image path without extension (for appending)
	
	// original image (all channels)
	orig_img = getTitle(); //save title of original image window to allow to open later for measures
		
	selectWindow(orig_img);
	getDimensions(ww, hh, cc, ss, ff);
	run("Clear Results");
	run("Set Measurements...", "area mean standard min integrated stack redirect=None decimal=3");
	nROI = roiManager("count");
	if(nROI == 0) {
		exit("No ROIs in ROI manager");
	} else {
		roiManager("Save", dir + File.separator + resultBaseName + "_RoiSet.zip");
	}
		
	for (i = 0; i < nROI; i++) { //loop over ROIs
		roiManager("select", i); //select current ROI
		
		for (j = 0; j < cc; j++) { //loop over channels
			roiManager("select", i); //select current ROI
			Stack.setChannel(j + 1); //set each channel 
			setCh_no = (j + 1); //channel number for naming file 	
			run("Line Width...", "line=" + lw); //increase line width
			run("Properties... ", "  width=" + lw); // increase line width
			//run("Area to Line"); //convert area ROI to line 
			// coordinates must be interpolated to get every pixel rather than "corners" of ROI
			run("Interpolate", "interval=1");
			getSelectionCoordinates(xpoints, ypoints);
			// make a list of z and t to match xpoints
			Stack.getPosition(channel, slice, frame);
			cpoints = newArray(xpoints.length);
			zpoints = newArray(xpoints.length);
			tpoints = newArray(xpoints.length);
			Array.fill(cpoints, channel);
			Array.fill(zpoints, slice);
			Array.fill(tpoints, frame);
			p = getProfile();
			// retrieve the bg value
			bg = bgSubtractValue(50);
			bgpoints = newArray(xpoints.length);
			Array.fill(bgpoints, bg);
			// retrieve the max value
			mx = imgMaxValue(100);
			mxpoints = newArray(xpoints.length);
			Array.fill(mxpoints, mx);
			// make a table of the data
			Table.create("Points");
			Table.setColumn("X", xpoints);
			Table.setColumn("Y", ypoints);
			Table.setColumn("C", cpoints);
			Table.setColumn("Z", zpoints);
			Table.setColumn("T", tpoints);
			Table.setColumn("Intensity", p);
			Table.setColumn("BG", bgpoints);
			Table.setColumn("MX", mxpoints);
			// the results will be saved with filename (no ext) + ROI number + channel
			// since ROI is unique to Z and T, we get all the data with just an i and j loop
			Table.save(fPath + "_" + i + "_" + setCh_no + "_data.txt");
			run("Clear Results");
			roiManager("deselect");
		}
	}
//	run("Close All");
}

function bgSubtractValue(size) {
	getDimensions(width, height, channels, slices, frames);
	bgvals = newArray(8);
	// starting positions for sq roi
	xpos = newArray(0, floor((width - size) / 2), width - size, 0, width - size, 0, floor((width - size) / 2), width - size);
	ypos = newArray(0, 0, 0, floor((height - size) / 2), floor((height - size) / 2), height - size, height - size, height - size);
	
	for(i = 0; i < bgvals.length; i ++) {
		makeRectangle(xpos[i], ypos[i], size, size);
		bgvals[i] = getValue("Mean");
	}
	Array.getStatistics(bgvals, min, max, mean, stdDev);
	
	return(min);
}

function imgMaxValue(size) {
	getDimensions(width, height, channels, slices, frames);
	
	makeRectangle(floor ((width - size) / 2), floor((height - size) / 2), size, size);
	val = getValue("Mean");
	
	return(val);
}

	
	
