"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var lcg_1 = require("../helpers/lcg");
var FFT = /** @class */ (function () {
    function FFT() {
    }
    FFT.log2 = function (n) {
        var log = 0;
        var k = 1;
        while (k < n) {
            k *= 2;
            log += 1;
        }
        if (n !== (1 << log))
            throw new Error("FFT: Data length is not a power of 2!: ".concat(n));
        return log;
    };
    FFT.bitreverse = function (data) {
        var _a, _b;
        var n = data.length / 2;
        var nm1 = n - 1;
        var i = 0, j = 0;
        while (i < nm1) {
            if (i < j) {
                var ii = i * 2;
                var jj = j * 2;
                _a = [data[jj], data[ii]], data[ii] = _a[0], data[jj] = _a[1];
                _b = [data[jj + 1], data[ii + 1]], data[ii + 1] = _b[0], data[jj + 1] = _b[1];
            }
            var k = n >> 1;
            while (k <= j) {
                j -= k;
                k >>= 1;
            }
            j += k;
            i += 1;
        }
    };
    FFT.transform_internal = function (data, direction) {
        if (data.length === 0)
            return;
        var n = data.length / 2;
        if (n === 1)
            return;
        var logn = FFT.log2(n);
        FFT.bitreverse(data);
        for (var bit = 0; bit < logn; bit++) {
            var dual = 1 << bit;
            var w_real = 1.0;
            var w_imag = 0.0;
            var theta = 2.0 * direction * Math.PI / (2.0 * dual);
            var s = Math.sin(theta);
            var t = Math.sin(theta / 2.0);
            var s2 = 2.0 * t * t;
            // a = 0
            for (var b = 0; b < n; b += 2 * dual) {
                var i = 2 * b;
                var j = 2 * (b + dual);
                var wd_real = data[j];
                var wd_imag = data[j + 1];
                data[j] = data[i] - wd_real;
                data[j + 1] = data[i + 1] - wd_imag;
                data[i] += wd_real;
                data[i + 1] += wd_imag;
            }
            // a = 1 .. (dual-1)
            for (var a = 1; a < dual; a++) {
                var tmp_real = w_real - s * w_imag - s2 * w_real;
                var tmp_imag = w_imag + s * w_real - s2 * w_imag;
                w_real = tmp_real;
                w_imag = tmp_imag;
                for (var b = 0; b < n; b += 2 * dual) {
                    var i = 2 * (b + a);
                    var j = 2 * (b + a + dual);
                    var z1_real = data[j];
                    var z1_imag = data[j + 1];
                    var wd_real = w_real * z1_real - w_imag * z1_imag;
                    var wd_imag = w_real * z1_imag + w_imag * z1_real;
                    data[j] = data[i] - wd_real;
                    data[j + 1] = data[i + 1] - wd_imag;
                    data[i] += wd_real;
                    data[i + 1] += wd_imag;
                }
            }
        }
    };
    FFT.prototype.transform = function (data) {
        FFT.transform_internal(data, -1);
    };
    FFT.prototype.inverse = function (data) {
        FFT.transform_internal(data, 1);
        var n = data.length / 2;
        var norm = 1.0 / n;
        for (var d = 0; d < data.length; d++) {
            data[d] *= norm;
        }
    };
    FFT.prototype.test = function (data) {
        var nd = data.length;
        var copy = data.slice();
        this.transform(data);
        //console.log('After transform:', data);
        this.inverse(data);
        //console.log('After inverse:', data);
        var diff = 0.0;
        for (var i = 0; i < nd; i++) {
            var d = data[i] - copy[i];
            diff += d * d;
        }
        return Math.sqrt(diff / nd);
    };
    FFT.prototype.make_random = function (n) {
        // Interleaved real/imag, like Rust
        var rand = new lcg_1.LCG(12345, 1345, 16645, 1013904);
        var nd = 2 * n;
        var data = [];
        for (var i = 0; i < nd; i++) {
            data.push(rand.nextDouble());
        }
        return data;
    };
    return FFT;
}());
function main() {
    var _a;
    var n = parseInt((_a = process.argv[2]) !== null && _a !== void 0 ? _a : "16", 10);
    var fft = new FFT();
    var data = fft.make_random(n);
    console.log("Typescript Specialized number FFT, n=" + n);
    for (var i = 0; i < 10; i++) {
        FFT.transform_internal(data, -1);
        FFT.transform_internal(data, 1);
        console.log("Loop ".concat(i, " done"));
    }
    // print array
    //console.log(data);
    //console.log(`n=${n} => RMS Error=${rms}`);
}
if (require.main === module) {
    main();
}
