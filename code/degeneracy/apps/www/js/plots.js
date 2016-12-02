var outputBinding = new Shiny.OutputBinding();

$.extend(outputBinding, {
  find: function(scope) {
    return $(scope).find('.d3graph');
  },
  renderValue: function(el, data) {  
    wrapper(el, data);
  },
});
Shiny.outputBindings.register(outputBinding);


function wrapper(el, data) {  
  var w = 500,
      h = 200,
      yaw = 0.5,
      pitch = 0.5, 
      drag = false;
  
  var plot_dat = null,
      surf_plot = null;

  if(data != null) {
    plot_dat = data.values;
    surf_plot = data.surf_plot;
  }

  
  d3.select(el).selectAll(".plotHolder").remove();

  if(surf_plot) {
  var selection = d3.select(el).selectAll("svg")
      .data(d3.keys(plot_dat));
  
  var enterSelection = selection.enter()
    .append("div")
    .attr("class", "plotHolder");
      
  var h3 = enterSelection
    .append("h5")
    .text(function(d) {return d;});
    
  var svg = enterSelection  
    .append("svg")
    .attr("width", w)
    .attr("height", h);
    
  selection.exit().remove();
  
  var group = svg.append("g");
  
  var md = group
    .data(d3.values(plot_dat));

    var surf = [] 
    md.each(function(d) {  
      surf.push(
        d3.select(this)
          .surface3D(w, h)
          .surfaceHeight(function(d){ 
              return d;
          })
          .surfaceColor(function(d){
              var c = d3.hsl((10*d + 100), 0.6, 0.5).rgb();
              return "rgb("+parseInt(c.r)+","+parseInt(c.g)+","+parseInt(c.b)+")";
          })
        );
    });
  
    d3.select(el).on("mousedown",function() {
      drag = [d3.mouse(this), yaw, pitch];
    })
    .on("mouseup",function() {
      drag = false;
    })
    .on("mousemove",function() {
      if(drag){            
        var mouse = d3.mouse(this);
        yaw = drag[1] - (mouse[0] - drag[0][0])/50;
        pitch = drag[2] + (mouse[1] - drag[0][1])/50;
        pitch = Math.max(-Math.PI/2, Math.min(Math.PI/2, pitch));
  
        for(i = 0; i < surf.length; i++) {
          surf[i].turntable(yaw, pitch);
        }      
      }
    });
  } 
  else if(!surf_plot) {
      
  var margin = {top: 30, right: 30, bottom: 40, left: 10}
  
  var x = d3.scale.linear()
      .range([0, w]);
  
  var y = d3.scale.linear()
      .range([h, 0]);
  
  var color = d3.scale.category10();
  
  var xAxis = d3.svg.axis()
      .scale(x)
      .orient("bottom");
  
  var yAxis = d3.svg.axis()
      .scale(y)
      .orient("left")
      .ticks(5);
  
var line = d3.svg.line()
    .interpolate("basis")
    .x(function(d) { return x(d.x); })
    .y(function(d) { return y(d.y); });
  
  var svg = d3.select(el).append("svg")
      .attr("width", w + margin.left + margin.right)
      .attr("height", h + margin.top + margin.bottom)
      .attr("class", "plotHolder")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  if(plot_dat != null){

    x.domain([
      d3.min(d3.values(plot_dat), function(c) { return d3.min(c.x) }),
      d3.max(d3.values(plot_dat), function(c) { return d3.max(c.x) })
    ]);

    y.domain([
      d3.min(d3.values(plot_dat), function(c) { return d3.min(c.value) }),
      d3.max(d3.values(plot_dat), function(c) { return d3.max(c.value) })
    ]);

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + h + ")")
        .call(xAxis)
      .append("text")
        .attr("transform", "translate(" + (w / 2) + " , " + parseFloat(margin.top + 0.5*margin.bottom) + ")")
        .style("text-anchor", "middle")
        .text(d3.keys(plot_dat)[0].split("_")[1]);
    
    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")
        .text("Expected Value");
    
    var lines = svg.selectAll(".vars")
      .data(d3.values(plot_dat))
      .enter()
      .append("g")
        .attr("class", "vars");
    

    lines.append("path")
      .attr("class", "line")
      .attr("d", function(d) {
        var tmp = []
        for(i =0 ; i < d.x.length; i++) tmp.push({x: d.x[i], y: d.value[i]});
        return line(tmp); 
      })
      .style("stroke", function(d, i) { return color(i); });
    
    var  legendSpace = w/d3.keys(plot_dat).length;
    
    lines.append("text")                            
          .attr("x", function(d,i){ return legendSpace/2+i*legendSpace })
          .attr("y", h + margin.top + margin.bottom) 
          .attr("class", "legend")   
          .style("fill", function(d, i) { return color(i) })
          .text(function(d, i) { return d3.keys(plot_dat)[i]; });  

    }
  }

}

  
