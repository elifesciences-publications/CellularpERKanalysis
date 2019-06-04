// Caroline Wee Jun 4 2019
// To extract GFP-positive cells from GFP channel
// Also extracts other cells using tERK channel
// First, open an image and crop/make substacks to isolate region of interest
// Use CropandSubstack_manualCW_v1.ijm
// Then, SAVE ALL CROPPED IMAGES TO A SEPARATE FOLDER
// Need to download latest FIJI version for Z-information to be accurate
// Do not click on images during analysis (if not in batch mode)!
// May need to be heavily optimized for your own images.
// Set Batchmode to "false" in order to visualize the process

// PARAMETERS THAT NEED TO BE MANUALLY SET
GFPchannel = 1;// need to change depending on image
pERKchannel = 2;
tERKchannel = 3; 
condition = "control" // treatment that fish is exposed to e.g. "shock" 

// PARAMETERS THAT NEED TO BE OPTIMIZED FOR YOUR IMAGES
ovalsize = 15; // sets number of pixels of circle diameter - need to change depending on image magnification)
minparticlesize = 1; // for tERK channel - sets threshold particle size that will be considered to be a "cell"
maxparticlesize = 25;
minparticlesize_GFP = 5; //same but for GFP channel
maxparticlesize_GFP = 30;
noisetolerance = 20; // increase if image is grainy

// PARAMETERS BELOW ARE IMPORTANT FOR PROPER THRESHOLDING
substack = 1; // yes = 1 -- basically whether there is GFP fluorescence in every image 
poortERK = 0; //if tERK quality is poor, additional processing will be implemented for cell segmentation

// Start of analysis
BatchMode = true; // specify "false" to verify it works - you will be able to visualize detected cells and procedure
setBatchMode(BatchMode); 

source_dir = getDirectory("Source Directory");
target_dir = getDirectory("Target Directory");

// make new directory
File.makeDirectory(target_dir + condition);

list = getFileList(source_dir);
list = Array.sort(list);

for (j=0; j<list.length; j++) {
    Image = source_dir + list[j];
    open(Image);
    name = getTitle();
    rename("Image");

    run("Clear Results");

    // This is to make a duplicate image that can then be stored for later use
    //while the tERK channel is processed for cell segmentation
    run("Duplicate...", "duplicate");
    rename("Imagefull");
    run("Duplicate...", "duplicate");
    rename("Imagefull2");
    selectWindow("Image");

    // Processing of tERK channel *need to close other two
    run("Split Channels");

    // close pERK and GFP channel
	close("C" + pERKchannel + "-Image");
  
    // Analysis of tERK
    selectWindow("C" + tERKchannel + "-Image");

    // Image processing for cell segmentation. These could be tweaked.
    run("Invert", "stack");
    //run("Macro...", "code = [v=" + maximum + "-v] stack");
    run("Subtract Background...", "rolling = 100 stack");
    run("Smooth", "stack");
    run("Gaussian Blur...", "sigma=2 stack");
    run("Enhance Contrast", "staurated=0.4");

    // Finding maxima for entire stack
    selectWindow("C" + tERKchannel + "-Image");
    input = getImageID();
    n = nSlices();
    k = 1;
        
    for (i=1; i<=n; i+=k) {
        showProgress(i, n);
        selectImage(input);
        setSlice(i);
        // May need to change noise threshold
        run("Find Maxima...", "noise=" + noisetolerance + " output=[Maxima Within Tolerance] exclude");
            
        if (i==1) 
            output = getImageID();
        else {
            run("Select All");
            run("Copy");
            close();
            selectImage(output);
            run("Add Slice");
            run("Paste");
         }
    }
  	
    run("Select None");

    // Run Analyze Particles to only extract particles of certain size and circularity (set below). I
    selectWindow("C" + tERKchannel + "-Image(1) Maxima");

    // Additional processing if tERK staining is poor
    if (poortERK ==1){
        run("Open", "stack");
        run("Erode", "stack");
        run("Watershed", "stack");
    }

 	// Here is where you find all particles of the correct size. Might need to modify
	run("Analyze Particles...", "size=" + minparticlesize + "-" + maxparticlesize + " circularity = 0.05-1.00 show=Outlines display stack");

    // Now to measure intensities of all extracted particles
    selectWindow("Imagefull");
    run("Set Measurements...", "mean center stack nan redirect=None decimal=3");

    npoints = nResults;
    xvalue = newArray(npoints);
    yvalue = newArray(npoints);
    zvalue = newArray(npoints);

    // to get the x y z coordinates for all points selected
 
    for (i=0; i<npoints; i+=k) {
        x = getResult("XM",i);
        y = getResult("YM", i);
        z = getResult("Slice", i);
        xvalue[i]= x;
        yvalue[i]= y;
        zvalue[i] = z;
	}

    run("Clear Results"); //now is when I make the actual measurements 
    run("Input/Output...", "jpeg=85 gif=-1 file =.xls");

    // to draw ovals over each point, and measure the intensity values
    for (i=0; i<npoints; i+=k) {
    	Stack.setChannel(1);
    	Stack.setSlice(zvalue[i]) 
    	toUnscaled(xvalue[i], yvalue[i]);
    	makeOval(xvalue[i]-ovalsize/2, yvalue[i]-ovalsize/2, ovalsize, ovalsize);
    	run("Measure");
    	fill();
    	setColor(0);

        // Measure pERK and tERK values
		Stack.setChannel(pERKchannel);
        run("Measure");
        Stack.setChannel(tERKchannel);
        run("Measure");
        fill();
        run("Select None");
        Stack.setChannel(1);
	}

	saveAs("Results",  target_dir + condition + '/' + name + "-allcells.xls");

	//To close open windows 
	run("Clear Results"); 
	if (BatchMode == false){
	selectWindow("C" + tERKchannel + "-Image(1) Maxima");
	run("Close");
	selectWindow("Drawing of " + "C" + tERKchannel + "-Image(1) Maxima");
	run("Close");
	selectWindow("Imagefull");
	run("Close");
	selectWindow("C" + tERKchannel + "-Image");
	run("Close");
	}

	//now for extracting cells using GFP channel

	//Image processing for cell segmentation
	selectWindow("C" + GFPchannel + "-Image");
	run("8-bit");

	// "calculate threshold for each image" necessary if there is signal in every image
	// that is, a substack is being used. 
	if (substack==1) {
	run("Make Binary", "method=Yen background=Dark calculate");
	} else {
	run("Make Binary", "method=Yen background=Dark");	
	}

	run("Open", "stack");
	run("Erode", "stack");
	run("Watershed", "stack");
	run("Analyze Particles...", "size=" + minparticlesize_GFP + "-" + maxparticlesize_GFP + " circularity = 0.05-1.00 show=Outlines display stack");

    // Now to measure intensities of all extracted particles
    selectWindow("Imagefull2");

    npoints = nResults;
    xvalue = newArray(npoints);
    yvalue = newArray(npoints);
    zvalue = newArray(npoints);

    // to get the x y z coordinates for all points selected
    for (i=0; i<npoints; i+=k) {
        x = getResult("XM",i);
        y = getResult("YM", i);
        z = getResult("Slice", i);
        xvalue[i]= x;
        yvalue[i]= y;
        zvalue[i] = z;
    }

    run("Clear Results"); // now is when I make the actual measurements 
    run("Input/Output...", "jpeg=85 gif=-1 file =.xls");

    // to draw ovals over each point, and measure the intensity values
    for (i=0; i<npoints; i+=k) {
		Stack.setChannel(GFPchannel); //GFP channel first
		Stack.setSlice(zvalue[i]) 
		toUnscaled(xvalue[i], yvalue[i]);
		makeOval(xvalue[i]-ovalsize/2, yvalue[i]-ovalsize/2, ovalsize, ovalsize);
		run("Measure");
		fill();
		setColor(0);

        // Measure pERK and tERK values
		Stack.setChannel(pERKchannel)
        run("Measure");
        Stack.setChannel(tERKchannel)
        run("Measure");
        fill();
        run("Select None");
        Stack.setChannel(1);
    }

    run("Close All");
    saveAs("Results",  target_dir + condition + '/' + name + "-GFP.xls");
}

setBatchMode(false);

