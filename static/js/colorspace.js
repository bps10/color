
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

    d3.select("#stimParam").on("change", changeStim);
    d3.select("#fundamentals").on("change", changeFund);
    d3.select("#Lpeak").on("change", changeLpeak);
    d3.select("#Mpeak").on("change", changeMpeak);
    d3.select("#Speak").on("change", changeSpeak);

    function changeStim() {
        stim = this.value;
        // post and then redraw w/ new values
    }

    function changeFund() {
        fundamental = this.value;
    }

    function changeLpeak() {
        Lpeak = this.value;

        d3.select("#Lpeak_val").remove();

        d3.select("#Lpeak_text").append("span")
        .attr("id", "Lpeak_val")
        .style("color", "red")
        .text(Lpeak + " nm");
    };

    function changeMpeak() {
        Mpeak = this.value;
        
        d3.select("#Mpeak_val").remove();
        
        d3.select("#Mpeak_text").append("span")
        .attr("id", "Mpeak_val")
        .style("color", "green")
        .text(Mpeak + " nm");
    };

    function changeSpeak() {
        Speak = this.value;
        
        d3.select("#Speak_val").remove();
        
        d3.select("#Speak_text").append("span")
        .attr("id", "Speak_val")
        .style("color", "blue")
        .text(Speak + " nm");
    };


