
// dimensions //
var margin = {top: 70, right: 20, bottom: 30, left: 40},
    width = 800 - margin.left - margin.right,
    height = 450 - margin.top - margin.bottom,
    xaxisWidth = 440, yaxisWidth = 340;


// scales //
var x1 = d3.scale.linear()
    .range([0, xaxisWidth]);

var y1 = d3.scale.linear()
    .range([yaxisWidth, 0]);

// axes handles //	
var xAxis = d3.svg.axis()
    .scale(x1)
    .ticks(6)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y1)
    .ticks(5)
    .orient("left");


var line = d3.svg.line()
	.x(function(d) { return x1(d.x); })
	.y(function(d) { return y1(d.y); });


function dataFormatter(x,y) {
    var dat = []
    for (var i=0; i<x.length; i++) {
        dat.push({"x": x[i], "y": y[i]});
        };
    return dat
};

var x = [0.1,0.2,0.3,0.4,0.5], 
	y = [0.5,0.4,0.3,0.2,0.1];
var data = dataFormatter(x, y);

var svg = d3.select("body").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	x1.domain([-0.3, 1.2]);
	y1.domain([-0.2, 1.2]);


    // figure 1 lines //

    svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(300," + (yaxisWidth - 60) + ")")
    .call(xAxis)
    .append("text")
    .attr("y", 40)
    .attr("x", 200)
    .style("text-anchor", "beginning")
    .text("r");


    svg.append("g")
    .attr("class", "y axis")
    .attr("transform", "translate(388, -11)")
    .call(yAxis)
    .append("text")
    .attr("transform", "rotate(-90)")
    .attr("y", -40)
    .attr("x", -150)
    .text("g");

    svg.append("path")
    .datum(data)
    .attr("transform", "translate(300, -11)")
    .attr("class", "line")
    .attr("d", line);

	///////////////////////
	// Input control box //
	///////////////////////
    var fundamentals = ["Neitz", "Stockman"],
        primaries = ["Wright", "CIE 1932", "Stiles and Burch"];

	d3.select("body").append("span")
		.attr("class", "inputholder");

    d3.select(".inputholder").append("span")
        .attr("id", "fundamentals")
        .text("fundamentals");
    d3.select("#fundamentals")
    .selectAll("li")
    .data(fundamentals)
    .enter()
    .append("li")
    //.attr("id", "fundamentals")
    //.attr("id", function(d) {return "fund_" + d;})
    .text(function(d) {return d;})
    .classed("selected", function(d) { return d === "Neitz";})
    .on("click", function(d) {
        xAxis = d;
        updateFund();
        });
    d3.select(".inputholder").append("br");

    d3.select(".inputholder").append("span")
        .attr("id", "primaries")
        .text("primaries");
    d3.select("#primaries")
    .selectAll("li")
    .data(primaries)
    .enter()
    .append("li")
    //.attr("id", "primaries")
    //.attr("id", function(d) {return "prim_" + d;})
    .text(function(d) {return d;})
    .classed("selected", function(d) { return d === "Wright";})
    .on("click", function(d) {
        xAxis = d;
        updatePrim();
        });
    d3.select(".inputholder").append("br");

	// L cone peak
	createSlider("inputholder", 559, "L");
	
	// M cone peak
	createSlider("inputholder", 530, "M");

	// S cone peak
	createSlider("inputholder", 417, "S");
	
    d3.select("#Lpeak").on("change", changeLpeak);
    d3.select("#Mpeak").on("change", changeMpeak);
    d3.select("#Speak").on("change", changeSpeak);

    function updatePrim() {
        console.log(xAxis);
        d3.select("#primaries")
        .selectAll("li")
        .classed("selected", function(d) {
                 return d === xAxis;
                 });}

    function updateFund() {
        console.log(xAxis);
        d3.select("#fundamentals")
        .selectAll("li")
        .classed("selected", function(d) {
                 return d === xAxis;
                 });}

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
		if (cone === "S") {color = "blue";}
		if (cone === "M") {color = "green";}
		if (cone === "L") {color = "red";}
		
        d3.select("#" + cone + "peak_val").remove();
        
        d3.select("#" + cone + "peak_text").append("span")
        .attr("id", cone + "peak_val")
        .style("color", color)
        .text(value + " nm");	
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
