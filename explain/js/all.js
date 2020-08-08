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