"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.LCG = void 0;
// Linear Congruential Generator (LCG)
var LCG = /** @class */ (function () {
    function LCG(seed, m, a, c) {
        this.last_num = seed;
        this.m = m;
        this.a = a;
        this.c = c;
    }
    LCG.prototype.nextDouble = function () {
        return this.nextInt() / this.m;
    };
    LCG.prototype.nextInt = function () {
        this.last_num = (this.a * this.last_num + this.c) % this.m;
        return this.last_num;
    };
    return LCG;
}());
exports.LCG = LCG;
