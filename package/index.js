var fs = require('fs')
var childProcess = require('child_process')

exports.handler = (event, context, callback) => {
    var file = null;
    if (event != null && event['z'] > 1) {
        console.log('A')
        file = fs.readFileSync('image1.jpg')
    }
    else if (event != null & event['z'] <= 1) {
        console.log('B')
        file = fs.readFileSync('image2.jpg')
    }
    if (fs.existsSync('index.js')) {
        console.log('index.js exists')
    }
    if (!fs.existsSync('hello.txt')) {
        console.log('hello.txt does not exist')
    }
    var binaryData = (new Buffer(file))
    var base64Data = binaryData.toString('base64')
    callback(null, base64Data);
};

// xxx = {
//     'z': '2',
//     'x': '1',
//     'y': '1'
// }

// function callback(a, b) {
//     var ls = childProcess.spawnSync('ps', ['-uxa'])
//     console.log(xxx['z'] + ' ' + xxx['x'] + ' ' + xxx['y'])
//     console.log(ls.stdout.toString())
// }

// exports.handler(xxx, null, callback)
