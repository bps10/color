

// dimensions //
var margin = {top: 70, right: 20, bottom: 30, left: 40},
    width = 850 - margin.left - margin.right,
    height = 450 - margin.top - margin.bottom,
    xaxisWidth = 490, yaxisWidth = 340;


// scales //
var x = d3.scale.linear()
    .range([0, xaxisWidth]);

var y = d3.scale.linear()
    .range([yaxisWidth, 0]);

// axes handles //	
var xAxis = d3.svg.axis()
    .scale(x)
    .ticks(6)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
    .ticks(5)
    .orient("left");


var line = d3.svg.line()
	.x(function(d) { return x(d.x); })
	.y(function(d) { return y(d.y); });


function dataFormatter(x,y) {
    var dat = []
    for (var i=0; i<x.length; i++) {
        dat.push({"x": x[i], "y": y[i]});
        };
    return dat
};

var L_lambdaMax = 559,
	M_lambdaMax = 530,
	S_lambdaMax = 420,
	opticalDensity = 0.2,
	output = "anti-log",
	startWavelength = 390,
	endWavelength = 750,
	step = 1;

var Ldata = dataFormatter(
				range(startWavelength, endWavelength, step),
				spectsens(L_lambdaMax, opticalDensity, output,
                startWavelength, endWavelength, step));

var Mdata = dataFormatter(
				range(startWavelength, endWavelength, step),
				spectsens(M_lambdaMax, opticalDensity, output,
                startWavelength, endWavelength, step));
				
var Sdata = dataFormatter(
				range(startWavelength, endWavelength, step),
				spectsens(S_lambdaMax, opticalDensity, output,
                startWavelength, endWavelength, step));
				
var svg = d3.select("body").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	x.domain(d3.extent(Ldata, function(d) { return d.x; }));
	y.domain(d3.extent(Ldata, function(d) { return d.y; }));


    // figure 1 lines //

    svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(300," + (yaxisWidth - 20) + ")")
    .call(xAxis)
    .append("text")
    .attr("y", 40)
    .attr("x", 200)
    .style("text-anchor", "beginning")
    .text("wavelength (nm)");


    svg.append("g")
    .attr("class", "y axis")
    .attr("transform", "translate(285, -40)")
    .call(yAxis)
    .append("text")
    .attr("transform", "rotate(-90)")
    .attr("y", 15)
    .attr("x", -70)
    .text("sensitivity");

    svg.append("path")
    .datum(Ldata)
    .attr("transform", "translate(300, -40)")
    .attr("class", "Lline")
    .attr("d", line);

    svg.append("path")
    .datum(Mdata)
    .attr("transform", "translate(300, -40)")
    .attr("class", "Mline")
    .attr("d", line);

    svg.append("path")
    .datum(Sdata)
    .attr("transform", "translate(300, -40)")
    .attr("class", "Sline")
    .attr("d", line);
	
	///////////////////////
	// Input control box //
	///////////////////////

	d3.select("body").append("span")
		.attr("class", "inputholder");
	// L cone peak
	createSlider("inputholder", 559, "L");
	
	// M cone peak
	createSlider("inputholder", 530, "M");

	// S cone peak
	createSlider("inputholder", 417, "S");
	
    d3.select("#Lpeak").on("change", changeLpeak);
    d3.select("#Mpeak").on("change", changeMpeak);
    d3.select("#Speak").on("change", changeSpeak);

    function changeLpeak() {
        Lpeak = this.value;
        updateSlider(Lpeak, "L");
    }

    function changeMpeak() {
        Mpeak = this.value;
        updateSlider(Mpeak, "M");
    }

    function changeSpeak() {
        Speak = this.value;       
        updateSlider(Speak, "S");
    }

	function updateSlider(value, coneType) {
		cone = coneType.toUpperCase();
		if (cone === "S") {color = "blue"; data = Sdata;}
		if (cone === "M") {color = "green"; data = Mdata;}
		if (cone === "L") {color = "red"; data = Ldata;}
		
        d3.select("#" + cone + "peak_val").remove();
		d3.select("." + cone + "line").remove();
        
        d3.select("#" + cone + "peak_text").append("span")
        .attr("id", cone + "peak_val")
        .style("color", color)
        .text(value + " nm");	
		
		data = dataFormatter(
						range(startWavelength, endWavelength, step),
						spectsens(value, opticalDensity, output,
						startWavelength, endWavelength, step));	
						
		svg.append("path")
		.datum(data)
		.attr("transform", "translate(300, -40)")
		.attr("class", cone + "line")
		.attr("d", line);
	}

	function createSlider(holder_name, value, coneType) {
	
		cone = coneType.toUpperCase();
		if (cone === "S") {color = "blue";}
		if (cone === "M") {color = "green";}
		if (cone === "L") {color = "red";}
		
		// S cone peak
		d3.select("." + holder_name)
			.append("input")
			.attr("id", cone + "peak")
			.attr("type", "range")
			.attr("min", 400)
			.attr("max", 600)
			.attr("step", 1)
			.attr("value", value);

		d3.select("." + holder_name).append("span")
			.attr("id", cone + "peak_text")
			.style("color", color)
			.text(" " + cone + ": ")
			.append("span")
			.attr("id", cone + "peak_val")
			.text(value + " nm");
		d3.select(".inputholder").append("br");	
	}