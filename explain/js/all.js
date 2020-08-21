var numPoints = 64;

function createRandomPoints(numPoints) {
	return Array.from(Array(numPoints)).map(function(a) {
		return Math.random() - .5;
	})
}
var datapoints1 = createRandomPoints(numPoints);
var datapoints2 = createRandomPoints(numPoints);
var step = 0; // initialization

var chartData = [{
		order: 1,
		type: 'bar',
		label: 'noise 1',
		data: datapoints1,
		borderColor: window.chartColors.blue,
		backgroundColor: window.chartColors.blue,
		//backgroundColor: 'rgba(0, 0, 0, 0)',
		//fill: false,
		lineTension: 0
	}, {
		order: 1,
		type: 'bar',
		label: 'noise 2',
		data: datapoints2,
		borderColor: window.chartColors.yellow,
		backgroundColor: window.chartColors.yellow,
		//backgroundColor: 'rgba(0, 0, 0, 0)',
		//fill: false,
		lineTension: 0
	},
	{
		order: 0,
		type: 'scatter',
		label: 'zero-crossings',
		backgroundColor: window.chartColors.red,
		data: [],
		fill: false,
		showLine: false,
		xAxisID: 'x-axis-2'
	}
];

function computeZeroCrossings2D() {
	chartData2D[1].data = [];
	var values = chartData2D[0].data;
	for (var j = 0; j < numPoints; j++) {
		for (var i = 0; i < numPoints; i++) {
			var idx0 = j * numPoints + i;
			var idx1 = j * numPoints + ((i + 1) >= numPoints ? 0 : i + 1);
			var idx2 = j * numPoints + ((i - 1) < 0 ? numPoints - 1 : i - 1);
			var idx3 = ((j + 1) >= numPoints ? 0 : j + 1) * numPoints + i;
			var idx4 = ((j - 1) < 0 ? numPoints - 1 : j - 1) * numPoints + i;
			if ((Math.sign(values[idx0].value) > 0 && Math.sign(values[idx1].value) < 0) || (Math.sign(values[idx0].value) < 0 && Math.sign(values[idx1].value) > 0)) {
				chartData2D[1].data.push({
					x: i + .5,
					y: j
				});
			}
			if ((Math.sign(values[idx0].value) > 0 && Math.sign(values[idx3].value) < 0) || (Math.sign(values[idx0].value) < 0 && Math.sign(values[idx3].value) > 0)) {
				chartData2D[1].data.push({
					x: i,
					y: j + .5
				});
			}
		}
	}
}

function computeZeroCrossings() {
	// recompute the zero crossings
	chartData[2].data = [];
	chartData[0].data.map(function(a, i, ar) {
		var cur = a;
		var next = ar[i + 1] || ar[0];
		if ((Math.sign(cur) > 0 && Math.sign(next) < 0) || (Math.sign(cur) < 0 && Math.sign(next) > 0)) {
			chartData[2].data.push({
				x: i + .5,
				y: 0
			});
		}
	});

}

computeZeroCrossings();

var config = {
	data: {
		labels: Array.from(Array(numPoints)).map(function(a, i) {
			return '' + i;
		}),
		datasets: chartData
	},
	options: {
		responsive: true,
		title: {
			display: true,
			text: 'Low-pass filtered white noise'
		},
		tooltips: {
			mode: 'index'
		},
		scales: {
			xAxes: [{
					display: true,
					scaleLabel: {
						display: true
					}
				},
				{
					id: 'x-axis-2',
					type: 'linear',
					display: false,
					ticks: {
						beginAtZero: true,
						min: 0,
						max: numPoints - 1
					}
				}
			],
			yAxes: [{
				display: false,
				scaleLabel: {
					display: true,
					labelString: '[AU]'
				},
			}]
		}
	}
};

var chartData2D = [{
	order: 1,
	type: 'scatter',
	label: 'noise',
	backgroundColor: window.chartColors.blue,
	data: [],
	pointBackgroundColor: [],
	fill: false,
	showLine: false,
	xAxisID: 'x-axis-2'
}, {
	order: 0,
	type: 'scatter',
	label: 'zero-crossings',
	backgroundColor: window.chartColors.red,
	data: [],
	fill: false,
	pointRadius: 1,
	showLine: false,
	xAxisID: 'x-axis-2'
}];

chartData2D[0].data = Array.from(Array(numPoints * numPoints)).map(function(a, i) {
	var y = Math.floor(i / numPoints);
	var x = i - (y * 64);
	var value = Math.random() - .5;
	var colormap = colorbrewer.Spectral[11];

	chartData2D[0].pointBackgroundColor.push(Math.random() - .5);
	return {
		value: value,
		x: x,
		y: y
	};
});

function values2Color() {
	var mm = Math.max(...chartData2D[0].data.map(function(a) {
		return Math.abs(a.value);
	}));
	chartData2D[0].pointBackgroundColor = Array.from(chartData2D[0].data).map(function(a, i) {
		var colormap = colorbrewer.Spectral[11];
		var value = ((a.value + mm) / (2 * mm)) * (colormap.length);
		return colormap[value.toFixed(0)];
	});
}
values2Color();

var config2D = {
	data: {
		datasets: chartData2D
	},
	options: {
		responsive: true,
		title: {
			display: true,
			text: 'Low-pass filtered white noise'
		},
		tooltips: {
			mode: 'index'
		},
		scales: {
			xAxes: [{
					display: true,
					scaleLabel: {
						display: true
					}
				},
				{
					id: 'x-axis-2',
					type: 'linear',
					display: false,
					ticks: {
						beginAtZero: true,
						min: 0,
						max: numPoints - 1
					}
				}
			],
			yAxes: [{
				display: false,
				scaleLabel: {
					display: true,
					labelString: '[AU]'
				},
			}]
		}
	}
}

var container, scene, camera, renderer, controls, stats;
var keyboard = new THREEx.KeyboardState();
var clock = new THREE.Clock();

// custom global variables
var points = [];
var values = [];
var points2 = [];
var values2 = [];
var size3d = 32;;

init3d();
animate();

function computeZeroCrossings3D(values, values2) {
    // do a very cheap computation not on the vertices but on the voxel given a threshold, return points
    var points = [];
    var axisMin = -10;
    var axisMax =  10;
    var axisRange = axisMax - axisMin;
    for (var k = 0; k < size3d; k++) {
	for (var j = 0; j < size3d; j++) {
	    for (var i = 0; i < size3d; i++) {
		var idx000 = i + (j * size3d) + (k * (size3d * size3d));
		if (Math.abs(values[idx000]) < 0.0005 && Math.abs(values2[idx000]) < 0.0005) {
		    var x = axisMin + axisRange * i / (size3d - 1);
		    var y = axisMin + axisRange * j / (size3d - 1);
		    var z = axisMin + axisRange * k / (size3d - 1);
		    
		    points.push( new THREE.Vector3(x, y, z));
		}
	    }
	}
    }    

    // add them to the scene
    var geometry = new THREE.CubeGeometry( .4, .4, .4 );
    var material = new THREE.MeshLambertMaterial( { color: 0xff5577 } );
    count = points.length; 
    var box = new THREE.InstancedMesh( geometry, material, count );
    box.name = "sphere";
    var transform = new THREE.Object3D();
    for (var i = 0; i < points.length; i++) {
	transform.position.set(points[i].x, points[i].y, points[i].z);
	transform.updateMatrix();
	box.setMatrixAt(i, transform.matrix);
	scene.add(box);
    }
    var oldBox = scene.getObjectByName("sphere");
    if (oldBox) {
	oldBox.geometry.dispose();
	oldBox.material.dispose();
	scene.remove(oldBox);
    }
    
    scene.add(box);
    return points;    
}

function filter3d( values, name, color) {
    // work on a copy
    var cV = [...values]; // values is now a local variable as argument to this function
    // average in 3 dimensions the values array based on the points array locations
    for (var k = 0; k < size3d; k++) {
	for (var j = 0; j < size3d; j++) {
	    for (var i = 0; i < size3d; i++) {
		var idx000 = i + (j * size3d) + (k * (size3d * size3d));
		ii = i+1; jj = j; kk = k;
		if (ii >= size3d)
		    ii = 0;
		if (ii < 0)
		    ii = size3d-1;
		if (jj >= size3d)
		    jj = 0;
		if (jj < 0)
		    jj = size3d-1;
		if (kk >= size3d)
		    kk = 0;
		if (kk < 0)
		    kk = size3d-1;
		var idx001 = ii + (jj * size3d) + (kk * (size3d * size3d));

		ii = i; jj = j+1; kk = k;
		if (ii >= size3d)
		    ii = 0;
		if (ii < 0)
		    ii = size3d-1;
		if (jj >= size3d)
		    jj = 0;
		if (jj < 0)
		    jj = size3d-1;
		if (kk >= size3d)
		    kk = 0;
		if (kk < 0)
		    kk = size3d-1;
		var idx010 = ii + (jj * size3d) + (kk * (size3d * size3d));
		
		ii = i; jj = j; kk = k+1;
		if (ii >= size3d)
		    ii = 0;
		if (ii < 0)
		    ii = size3d-1;
		if (jj >= size3d)
		    jj = 0;
		if (jj < 0)
		    jj = size3d-1;
		if (kk >= size3d)
		    kk = 0;
		if (kk < 0)
		    kk = size3d-1;
		var idx100 = ii + (jj * size3d) + (kk * (size3d * size3d));
		
		ii = i-1; jj = j; kk = k;
		if (ii >= size3d)
		    ii = 0;
		if (ii < 0)
		    ii = size3d-1;
		if (jj >= size3d)
		    jj = 0;
		if (jj < 0)
		    jj = size3d-1;
		if (kk >= size3d)
		    kk = 0;
		if (kk < 0)
		    kk = size3d-1;
		var idx00_1 = ii + (jj * size3d) + (kk * (size3d * size3d));

		ii = i; jj = j-1; kk = k;
		if (ii >= size3d)
		    ii = 0;
		if (ii < 0)
		    ii = size3d-1;
		if (jj >= size3d)
		    jj = 0;
		if (jj < 0)
		    jj = size3d-1;
		if (kk >= size3d)
		    kk = 0;
		if (kk < 0)
		    kk = size3d-1;
		var idx0_10 = ii + (jj * size3d) + (kk * (size3d * size3d));
		
		ii = i; jj = j; kk = k-1;
		if (ii >= size3d)
		    ii = 0;
		if (ii < 0)
		    ii = size3d-1;
		if (jj >= size3d)
		    jj = 0;
		if (jj < 0)
		    jj = size3d-1;
		if (kk >= size3d)
		    kk = 0;
		if (kk < 0)
		    kk = size3d-1;
		var idx_100 = ii + (jj * size3d) + (kk * (size3d * size3d));

		// just the nearest neighbors on a torus averaged - will produce edge effects
		cV[idx000] = (values[idx000] + values[idx001] + values[idx00_1] + values[idx010] + values[idx0_10] + values[idx100] + values[idx_100]) / 7.0;
	    }
	}
    }
    values = cV;
    // after filtering create another IsoSurface and display
    createIsoSurface( values, name, color );
    return values;
}

function createIsoSurface( values, name, color ) {
    // Marching Cubes Algorithm
    
    var size2 = size3d * size3d;
    
    // Vertices may occur along edges of cube, when the values at the edge's endpoints
    //   straddle the isolevel value.
    // Actual position along edge weighted according to function values.
    var vlist = new Array(12);
    
    var geometry = new THREE.Geometry();
    var vertexIndex = 0;
    
    for (var z = 0; z < size3d - 1; z++)
	for (var y = 0; y < size3d - 1; y++)
	    for (var x = 0; x < size3d - 1; x++)
    {
	// index of base point, and also adjacent points on cube
	var p    = x + size3d * y + size2 * z,
	    px   = p   + 1,
	    py   = p   + size3d,
	    pxy  = py  + 1,
	    pz   = p   + size2,
	    pxz  = px  + size2,
	    pyz  = py  + size2,
	    pxyz = pxy + size2;
	
	// store scalar values corresponding to vertices
	var value0 = values[ p    ],
	    value1 = values[ px   ],
	    value2 = values[ py   ],
	    value3 = values[ pxy  ],
	    value4 = values[ pz   ],
	    value5 = values[ pxz  ],
	    value6 = values[ pyz  ],
	    value7 = values[ pxyz ];
	
	// place a "1" in bit positions corresponding to vertices whose
	//   isovalue is less than given constant.
	
	var isolevel = 0;
	
	var cubeindex = 0;
	if ( value0 < isolevel ) cubeindex |= 1;
	if ( value1 < isolevel ) cubeindex |= 2;
	if ( value2 < isolevel ) cubeindex |= 8;
	if ( value3 < isolevel ) cubeindex |= 4;
	if ( value4 < isolevel ) cubeindex |= 16;
	if ( value5 < isolevel ) cubeindex |= 32;
	if ( value6 < isolevel ) cubeindex |= 128;
	if ( value7 < isolevel ) cubeindex |= 64;
	
	// bits = 12 bit number, indicates which edges are crossed by the isosurface
	var bits = THREE.edgeTable[ cubeindex ];
	
	// if none are crossed, proceed to next iteration
	if ( bits === 0 ) continue;
	
	// check which edges are crossed, and estimate the point location
	//    using a weighted average of scalar values at edge endpoints.
	// store the vertex in an array for use later.
	var mu = 0.5;
	
	// bottom of the cube
	if ( bits & 1 )
	{
	    mu = ( isolevel - value0 ) / ( value1 - value0 );
	    vlist[0] = points[p].clone().lerp( points[px], mu );
	}
	if ( bits & 2 )
	{
	    mu = ( isolevel - value1 ) / ( value3 - value1 );
	    vlist[1] = points[px].clone().lerp( points[pxy], mu );
	}
	if ( bits & 4 )
	{
	    mu = ( isolevel - value2 ) / ( value3 - value2 );
	    vlist[2] = points[py].clone().lerp( points[pxy], mu );
	}
	if ( bits & 8 )
	{
	    mu = ( isolevel - value0 ) / ( value2 - value0 );
	    vlist[3] = points[p].clone().lerp( points[py], mu );
	}
	// top of the cube
	if ( bits & 16 )
	{
	    mu = ( isolevel - value4 ) / ( value5 - value4 );
	    vlist[4] = points[pz].clone().lerp( points[pxz], mu );
	}
	if ( bits & 32 )
	{
	    mu = ( isolevel - value5 ) / ( value7 - value5 );
	    vlist[5] = points[pxz].clone().lerp( points[pxyz], mu );
	}
	if ( bits & 64 )
	{
	    mu = ( isolevel - value6 ) / ( value7 - value6 );
	    vlist[6] = points[pyz].clone().lerp( points[pxyz], mu );
	}
	if ( bits & 128 )
	{
	    mu = ( isolevel - value4 ) / ( value6 - value4 );
	    vlist[7] = points[pz].clone().lerp( points[pyz], mu );
	}
	// vertical lines of the cube
	if ( bits & 256 )
	{
	    mu = ( isolevel - value0 ) / ( value4 - value0 );
	    vlist[8] = points[p].clone().lerp( points[pz], mu );
	}
	if ( bits & 512 )
	{
	    mu = ( isolevel - value1 ) / ( value5 - value1 );
	    vlist[9] = points[px].clone().lerp( points[pxz], mu );
	}
	if ( bits & 1024 )
	{
	    mu = ( isolevel - value3 ) / ( value7 - value3 );
	    vlist[10] = points[pxy].clone().lerp( points[pxyz], mu );
	}
	if ( bits & 2048 )
	{
	    mu = ( isolevel - value2 ) / ( value6 - value2 );
	    vlist[11] = points[py].clone().lerp( points[pyz], mu );
	}
	
	// construct triangles -- get correct vertices from triTable.
	var i = 0;
	cubeindex <<= 4;  // multiply by 16...
	// "Re-purpose cubeindex into an offset into triTable."
	//  since each row really isn't a row.
	
	// the while loop should run at most 5 times,
	//   since the 16th entry in each row is a -1.
	while ( THREE.triTable[ cubeindex + i ] != -1 )
	{
	    var index1 = THREE.triTable[cubeindex + i];
	    var index2 = THREE.triTable[cubeindex + i + 1];
	    var index3 = THREE.triTable[cubeindex + i + 2];
	    
	    geometry.vertices.push( vlist[index1].clone() );
	    geometry.vertices.push( vlist[index2].clone() );
	    geometry.vertices.push( vlist[index3].clone() );
	    var face = new THREE.Face3(vertexIndex, vertexIndex+1, vertexIndex+2);
	    geometry.faces.push( face );
	    
	    geometry.faceVertexUvs[ 0 ].push( [ new THREE.Vector2(0,0), new THREE.Vector2(0,1), new THREE.Vector2(1,1) ] );
	    
	    vertexIndex += 3;
	    i += 3;
	}
    }
    
    //geometry.computeCentroids();
    geometry.computeFaceNormals();
    geometry.computeVertexNormals();
    
    var colorMaterial =  new THREE.MeshLambertMaterial( {color: color, side: THREE.DoubleSide, transparent: false, opacity: 0.9} );
    var oldIsoSurface = scene.getObjectByName(name);
    if (oldIsoSurface) { // remove the old surface
	oldIsoSurface.geometry.dispose();
	oldIsoSurface.material.dispose();
	scene.remove(oldIsoSurface);
	if (typeof renderer.renderLists !== 'undefined')
	    renderer.renderLists.dispose();
	//animate();
    }
    // add a new surface
    var IsoSurface = new THREE.Mesh( geometry, colorMaterial );
    IsoSurface.name = name;
    scene.add(IsoSurface);
}

function init3d() {
    scene = new THREE.Scene();
    // CAMERA
    var SCREEN_WIDTH = jQuery('#ThreeJS').width(), SCREEN_HEIGHT = jQuery('#ThreeJS').height();
    var VIEW_ANGLE = 45, ASPECT = SCREEN_WIDTH / SCREEN_HEIGHT, NEAR = 0.1, FAR = 20000;
    camera = new THREE.PerspectiveCamera( VIEW_ANGLE, ASPECT, NEAR, FAR);
    scene.add(camera);
    camera.position.set(20,20,60);
    camera.lookAt(scene.position);	
    // RENDERER
    if ( Detector.webgl )
	renderer = new THREE.WebGLRenderer( {antialias:true} );
    else
	renderer = new THREE.CanvasRenderer();
    renderer.setClearColor("#FFFFFF");
    renderer.setSize(jQuery('#ThreeJS').width(), jQuery('#ThreeJS').height());
    container = document.getElementById( 'ThreeJS' );
    container.appendChild( renderer.domElement );
    // EVENTS
    THREEx.WindowResize(renderer, camera);
    //THREEx.FullScreen.bindKey({ charCode : 'm'.charCodeAt(0) });
    // CONTROLS
    controls = new THREE.OrbitControls( camera, renderer.domElement );
    // STATS
    //stats = new Stats();
    //stats.domElement.style.position = 'absolute';
    //stats.domElement.style.bottom = '0px';
    //stats.domElement.style.zIndex = 100;
    //container.appendChild( stats.domElement );
    // LIGHT
    var light = new THREE.PointLight(0xffffff);
    light.position.set(0,50,0);
    scene.add(light);
    var light2 = new THREE.PointLight(0xffffff);
    light2.position.set(0,-50,0);
    scene.add(light2);
    
    ////////////
    // CUSTOM //
    ////////////
    
    scene.add( new THREE.AxisHelper(100) );
    
    // number of cubes along a side
    //size = 32;
    
    var axisMin = -10;
    var axisMax =  10;
    var axisRange = axisMax - axisMin;
    
    // Generate a list of 3D points and values at those points
    for (var k = 0; k < size3d; k++)
	for (var j = 0; j < size3d; j++)
	    for (var i = 0; i < size3d; i++)
    {
	// actual values
	var x = axisMin + axisRange * i / (size3d - 1);
	var y = axisMin + axisRange * j / (size3d - 1);
	var z = axisMin + axisRange * k / (size3d - 1);
	points.push( new THREE.Vector3(x,y,z) );
	points2.push( new THREE.Vector3(x,y,z) );
	var value = Math.random() - .5;// x*x + y*y - z*z - 25;
	values.push( value );
	var value = Math.random() - .5;// x*x + y*y - z*z - 25;
	values2.push( value );
    }
    createIsoSurface(values, "IsoSurface", "#5577ff");
    createIsoSurface(values2, "IsoSurface2", "#ffff00");
}

function animate() {
    requestAnimationFrame( animate );
    render();		
    update();
}

function update() {
    controls.update();
    if (typeof stats !== 'undefined')
	stats.update();
}

function render()  {
    renderer.render( scene, camera );
}


window.onload = function() {
    var ctx = document.getElementById('canvas').getContext('2d');
    window.myLine = new Chart(ctx, config);
    var ctx2D = document.getElementById('canvas2D').getContext('2d');
    window.myLine2D = new Chart(ctx2D, config2D);
    setTimeout(function() {
	chartData[1].hidden = true;
    }, 100);
    setTimeout(function() {
	chartData2D[1].hidden = true;
    }, 100);
    
    // 3D case
    //init3d();
};

function smooth2D(values) {
	var smoothed = [];
	for (var j = 0; j < numPoints; j++) {
		for (var i = 0; i < numPoints; i++) {
			var idx0 = j * numPoints + i;
			var idx1 = j * numPoints + ((i + 1) >= numPoints ? 0 : i + 1);
			var idx2 = j * numPoints + ((i - 1) < 0 ? numPoints - 1 : i - 1);
			var idx3 = ((j + 1) >= numPoints ? 0 : j + 1) * numPoints + i;
			var idx4 = ((j - 1) < 0 ? numPoints - 1 : j - 1) * numPoints + i;
			var improved = Number(this.average([values[idx0].value, values[idx1].value, values[idx2].value, values[idx3].value, values[idx4].value]));
			smoothed.push({
				x: values[idx0].x,
				y: values[idx0].y,
				value: improved
			});
		}
	}
	return smoothed;
}

function smooth(values, alpha) {
	var weighted = average(values) * alpha;
	var smoothed = [];
	for (var i in values) {
		var curr = values[i];
		var prev = values[i - 1] || values[values.length - 1];
		var next = values[i + 1] || values[0];
		var improved = Number(this.average([.5555 * prev, .888 * curr, .5555 * next]));
		smoothed.push(improved);
	}
	return smoothed;
}

function average(data) {
	var sum = data.reduce(function(sum, value) {
		return sum + value;
	}, 0);
	var avg = sum / data.length;
	return avg;
}
step2D = 0;
document.getElementById('filterData2').addEventListener('click', function() {
	if (step2D < 5) {
		chartData2D[0].data = smooth2D(chartData2D[0].data);
	} else if (step2D < 17) {
		for (var i = 0; i < 7; i++) {
			chartData2D[0].data = smooth2D(chartData2D[0].data);
		}
	}
	values2Color();
	computeZeroCrossings2D();
	step2D++;
	window.myLine2D.update();
});
step3D = 0;
document.getElementById('filterData3').addEventListener('click', function() {
    if (step3D < 5) {
	values = filter3d(values, "IsoSurface", "#5577ff");
	values2 = filter3d(values2, "IsoSurface2", "#ffff22");
    } else if (step3D < 17) {
	for (var i = 0; i < 7; i++) {
	    values = filter3d(values, "IsoSurface", "#5577ff");
	    values2 = filter3d(values2, "IsoSurface2", "#ffff00");
	}
    } 
    computeZeroCrossings3D(values, values2);
    step3D++;
    //window.myLine2D.update();
});


document.getElementById('randomizeData').addEventListener('click', function() {
	if (step == 0) {
		// smooth with Gaussian
		chartData[0].data = smooth(chartData[0].data, 0.85);
		chartData[1].data = smooth(chartData[1].data, 0.85);
		// put data back to mean 0
		var mean0 = average(chartData[0].data);
		chartData[0].data = chartData[0].data.map(function(a) {
			return a - mean0;
		});
		var mean1 = average(chartData[1].data);
		chartData[1].data = chartData[1].data.map(function(a) {
			return a - mean1;
		});
	} else if (step < 7) {
		console.log("smooth the two sets of points independently from each other...");
		for (var i = 0; i < 7; i++) {
			// smooth with Gaussian
			chartData[0].data = smooth(chartData[0].data, 0.0);
			chartData[1].data = smooth(chartData[1].data, 0.0);
			// put data back to mean 0
			var mean0 = average(chartData[0].data);
			chartData[0].data = chartData[0].data.map(function(a) {
				return a - mean0;
			});
			var mean1 = average(chartData[1].data);
			chartData[1].data = chartData[1].data.map(function(a) {
				return a - mean1;
			});
		}
	}
	step = step + 1;

	computeZeroCrossings();

	/* 	for (var i = 0, len = datapoints.length; i < len; ++i) {
			datapoints[i] = Math.random() < 0.05 ? NaN : randomScalingFactor();
        } */
	window.myLine.update();
});
