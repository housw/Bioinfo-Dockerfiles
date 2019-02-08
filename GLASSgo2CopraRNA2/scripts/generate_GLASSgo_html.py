#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
import json
import textwrap
import argparse


class HTML:

    def __init__(self):

        self.treeData = None

        self.head = textwrap.dedent("""
                    <!DOCTYPE html>
                    <html lang="en">
                      <head>
                        <meta charset="utf-8">

                        <script src="https://cdn.rawgit.com/eligrey/canvas-toBlob.js/f1a01896135ab378aa5c0118eadd81da55e698d8/canvas-toBlob.js"> </script>
                        <script src="https://cdn.rawgit.com/eligrey/FileSaver.js/e9d941381475b5df8b7d7691013401e171014e89/FileSaver.min.js"></script>

                        <title>Tree Example</title>

                        <style>

                            html,body{
                                height:100%;
                                margin:0;
                                overflow:hidden;
                            }

                            svg{
                                width:100%;
                                height:100%;
                            }

                            input[type="file"]{
                                position:fixed;
                                top:1em;
                                left:1em;
                            }

                            .node {
                                cursor: pointer;
                            }

                            .node circle {
                              fill: #fff;
                              stroke: steelblue;
                              stroke-width: 3px;
                            }

                            .node text {
                              font: 12px sans-serif;
                            }

                            .link {
                              fill: none;
                              stroke: #ccc;
                              stroke-width: 2px;
                            }

                        </style>

                      </head>
                    """)

        self.body_first = textwrap.dedent("""

                      <body>

                        <!--<div>
                            <button id='saveButton'>Export my D3 visualization to PNG</button>
                        </div> -->

                        <!-- Dateiupload -->
                        <!-- <input id="get-tree-data" type='file' accept='application/json' /> -->

                        <div id="treeviz"
                            style="width:100%; height:100%;border:solid;border-width:1px;">
                        </div>

                        <!-- load the d3.js library -->
                        <script src="http://d3js.org/d3.v3.min.js"></script>

                        <script>

                            var treeData =
                    """)

        self.body_second = textwrap.dedent("""
                        var margin = {top: 20, right: 120, bottom: 20, left: 120},
                            width = 650 - margin.right - margin.left,
                            height = 500 - margin.top - margin.bottom;


                        var i = 0,
                            duration = 750,
                            rectW = 60,
                            rectH = 30;
                            //root;

                        var tree = d3.layout.tree()
                            .separation(function separation(a,b) {
                                return a.parent == b.parent ? 1 : 0.5;
                            })
                            .nodeSize([50, 0]);


                        var diagonal = d3.svg.diagonal()
                            .projection(function(d) {
                                return [d.y, d.x];
                            });


                        var zoom = d3.behavior.zoom()
                            .scaleExtent([0.1, 10])
                            .on("zoom", zoomed);

                        var drag = d3.behavior.drag()
                            .origin(function(d) { return d; })
                            .on("dragstart", dragstarted)
                            .on("drag", dragged)
                            .on("dragend", dragended);


                        var svg = d3.select("#treeviz").append("svg")
                            .attr("width", width + margin.right + margin.left)
                            .attr("height", height + margin.top + margin.bottom)
                            .call(zoom)
                            .append("g")
                            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

                        root = treeData[0];
                        root.x0 = height;
                        root.y0 = 0;
                        update(root);
                        d3.select(self.frameElement).style("height", "500px");

                        function update(source) {

                          // Compute the new tree layout.
                          var nodes = tree.nodes(root).reverse(),
                              links = tree.links(nodes);

                          // Normalize for fixed-depth.
                          nodes.forEach(function(d) { d.y = d.depth * 320; });

                          // Update the nodes…
                          var node = svg.selectAll("g.node")
                              .data(nodes, function(d) { return d.id || (d.id = ++i); });

                          // Enter any new nodes at the parent's previous position.
                          var nodeEnter = node.enter().append("g")
                              .attr("class", "node")
                              .attr("transform", function(d) { return "translate(" + source.y0 + "," + source.x0 + ")"; })
                              .on("click", click);

                          nodeEnter.append("circle")
                           .attr("r", function(d) { return d.value; })
                           .style("stroke", function(d) { return d.type; })
                           .style("fill", function(d) { return d.level; });

                          nodeEnter.append("text")
                              .attr("x", function(d) { return d.children || d._children ? -5 : -5; })
                              .attr("dy", "0.35em")
                              //.attr("text-anchor", function(d) { return d.children || d._children ? "end" : "start"; })
                              .text(function(d) { return d.name; })
                              .style("fill-opacity", 1e-6);

                          // Transition nodes to their new position.
                          var nodeUpdate = node.transition()
                              .duration(duration)
                              .attr("transform", function(d) { return "translate(" + d.y + "," + d.x + ")"; });

                          nodeUpdate.select("circle")
                            .attr("r", function(d) { return d.value; })
                            .style("fill", function(d) { return d.level; });

                          nodeUpdate.select("text")
                              .style("fill-opacity", 1);

                          // Transition exiting nodes to the parent's new position.
                          var nodeExit = node.exit().transition()
                              .duration(duration)
                              .attr("transform", function(d) { return "translate(" + source.y + "," + source.x + ")"; })
                              .remove();

                          nodeExit.select("circle")
                              .attr("r", 1e-6);

                          nodeExit.select("text")
                              .style("fill-opacity", 1e-6);

                          // Update the links
                          var link = svg.selectAll("path.link")
                              .data(links, function(d) { return d.target.id; });

                          // Enter any new links at the parent's previous position.
                          link.enter().insert("path", "g")
                              .attr("class", "link")
                              .style("stroke", function(d) { return d.target.level; })
                              .attr("d", function(d) {
                                var o = {x: source.x0, y: source.y0};
                                return diagonal({source: o, target: o});
                              });

                          // Transition links to their new position.
                          link.transition()
                              .duration(duration)
                              .attr("d", diagonal);

                          // Transition exiting nodes to the parent's new position.
                          link.exit().transition()
                              .duration(duration)
                              .attr("d", function(d) {
                                var o = {x: source.x, y: source.y};
                                return diagonal({source: o, target: o});
                              })
                              .remove();

                          // Stash the old positions for transition.
                          nodes.forEach(function(d) {
                            d.x0 = d.x;
                            d.y0 = d.y;
                          });
                        }

                        // Toggle children on click.
                        function click(d) {
                          if (d.children) {
                            d._children = d.children;
                            d.children = null;
                          } else {
                            d.children = d._children;
                            d._children = null;
                          }
                          update(d);
                        }


                        function zoomed() {
                            svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
                        }

                        function dragstarted(d) {
                            d3.event.sourceEvent.stopPropagation();
                            d3.select(this).classed("dragging", true);
                        }

                        function dragged(d) {
                            d3.select(this).attr("x", d.x = d3.event.x).attr("y", d.y = d3.event.y);
                        }

                        function dragended(d) {
                            d3.select(this).classed("dragging", false);
                        }

                    if(window.File&&window.FileReader&&window.FileList){
                        function handleFileSelect(event){
                            var file=event.target.files[0];
                            var reader=new FileReader();
                            reader.onload=function(){
                                run(JSON.parse(reader.result));}
                            reader.onerror=function(){
                                alert("Too bad - error on loading file: "+file.name);}
                            reader.readAsText(file);}
                        document.getElementById('get-tree-data').
                            addEventListener('change',handleFileSelect,false);}
                    else{alert('Dude - your browser is too old! Do something!');}

                    // Set-up the export button
                    d3.select('#saveButton').on('click', function(){
                        var svgcanvas = document.getElementsByTagName("svg");
                        var xml = createXML(svgcanvas);

                        var height = parseInt(svgcanvas[0].height);
                        var width = parseInt(svgcanvas[0].width);

                    });

                    </script>
                  </body>
                </html>
                    """)

    def load_json(self, json_file):
        """
        :param json_file: input json file generated by GLASSgo
        :return: None
        """
        with open(json_file) as json_data:
            d = json.load(json_data)
            self.treeData = d

    def write_html(self, output_file):
        """
        :param output_file: write html file
        :return: None
        """

        with open(output_file, 'w') as oh:
            oh.write(self.head)

        with open(output_file, 'a') as oh:
            oh.write(self.body_first)
            json.dump(self.treeData, oh, indent=2, sort_keys=False)
            oh.write(";")
            oh.write(self.body_second)


def main():

    # main parser
    parser = argparse.ArgumentParser(description="visualize GLASSgo json to html") 
    parser.add_argument("input_json_file", help="input json file")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory, default=./", default="./")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    if args.prefix:
        prefix = args.prefix
        out_file = os.path.join(args.out_folder, prefix+".html")
    else:
        basename = os.path.basename(args.input_json_file)
        prefix = os.path.splitext(basename)[0]
        out_file = os.path.join(args.out_folder, prefix+".html")


    # convert
    html = HTML()
    html.load_json(json_file=args.input_json_file)
    html.write_html(output_file = out_file)


if __name__ == '__main__':
    main()
