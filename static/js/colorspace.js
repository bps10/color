
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
	d3.select("body").append("span")
		.attr("class", "inputholder");

	d3.select(".inputholder").append("ul").attr("id","fundamentals")
		.style("font-weight","bold")
		.text("Fundamentals")
	d3.select("#fundamentals").append("li")
		.attr("value","fund_neitz")
		.style("font-weight","normal")
		.text("Neitz");
	d3.select("#fundamentals")
		.attr("value","fund_stockman")
		.append("li")
		.style("font-weight","normal")
		.text("Stockman");

	d3.select(".inputholder").append("ul")
		.style("font-weight","bold")
		.attr("id","primaries")
		.text("Primaries")
	d3.select("#primaries").append("li")
		.style("font-weight","normal")
		.text("Wright");
	d3.select("#primaries").append("li")
		.style("font-weight","normal")
		.text("Stiles and Burch");
	d3.select("#primaries").append("li")
		.style("font-weight","normal")
		.text("CIE 1932");
	
	// L cone peak
	createSlider("inputholder", 559, "L");
	
	// M cone peak
	createSlider("inputholder", 530, "M");

	// S cone peak
	createSlider("inputholder", 417, "S");
	
	d3.select("#fundamentals").on("click", changeFund);
	d3.select("#primaries").on("click", changePrim);
    //d3.select("#stimParam").on("change", changeStim);
    //d3.select("#fundamentals").on("change", changeFund);
    d3.select("#Lpeak").on("change", changeLpeak);
    d3.select("#Mpeak").on("change", changeMpeak);
    d3.select("#Speak").on("change", changeSpeak);

	
	
    function changePrim() {
        stim = this.value;
		console.log(stim);
        // post and then redraw w/ new values
		// change class = selected
    }

    function changeFund() {
        fundamental = this.value;
		console.log(fundamental);
	}

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
/*
  // Build menus
  d3.select('#x-axis-menu')
    .selectAll('li')
    .data(xAxisOptions)
    .enter()
    .append('li')
    .text(function(d) {return d;})
    .classed('selected', function(d) {
      return d === xAxis;
    })
    .on('click', function(d) {
      xAxis = d;
      updateChart();
      updateMenus();
    });
	
  function updateMenus() {
    d3.select('#x-axis-menu')
      .selectAll('li')
      .classed('selected', function(d) {
        return d === xAxis;
      });
    d3.select('#y-axis-menu')
      .selectAll('li')
      .classed('selected', function(d) {
        return d === yAxis;
    });
  }
  */