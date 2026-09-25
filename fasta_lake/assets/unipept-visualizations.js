//#region src/DataNode.ts
var e = class {
	constructor(e, t, n = [], r, i, a = {}) {
		this.id = e, this.name = t, this.children = n, this.count = r, this.selfCount = i, this.extra = a;
	}
};
//#endregion
//#region node_modules/d3-array/src/ascending.js
function t(e, t) {
	return e == null || t == null ? NaN : e < t ? -1 : e > t ? 1 : e >= t ? 0 : NaN;
}
//#endregion
//#region node_modules/d3-array/src/descending.js
function n(e, t) {
	return e == null || t == null ? NaN : t < e ? -1 : t > e ? 1 : t >= e ? 0 : NaN;
}
//#endregion
//#region node_modules/d3-array/src/bisector.js
function r(e) {
	let r, a, o;
	e.length === 2 ? (r = e === t || e === n ? e : i, a = e, o = e) : (r = t, a = (n, r) => t(e(n), r), o = (t, n) => e(t) - n);
	function s(e, t, n = 0, i = e.length) {
		if (n < i) {
			if (r(t, t) !== 0) return i;
			do {
				let r = n + i >>> 1;
				a(e[r], t) < 0 ? n = r + 1 : i = r;
			} while (n < i);
		}
		return n;
	}
	function c(e, t, n = 0, i = e.length) {
		if (n < i) {
			if (r(t, t) !== 0) return i;
			do {
				let r = n + i >>> 1;
				a(e[r], t) <= 0 ? n = r + 1 : i = r;
			} while (n < i);
		}
		return n;
	}
	function l(e, t, n = 0, r = e.length) {
		let i = s(e, t, n, r - 1);
		return i > n && o(e[i - 1], t) > -o(e[i], t) ? i - 1 : i;
	}
	return {
		left: s,
		center: l,
		right: c
	};
}
function i() {
	return 0;
}
//#endregion
//#region node_modules/d3-array/src/number.js
function a(e) {
	return e === null ? NaN : +e;
}
//#endregion
//#region node_modules/d3-array/src/bisect.js
var o = r(t), s = o.right;
o.left, r(a).center;
//#endregion
//#region node_modules/internmap/src/index.js
var c = class extends Map {
	constructor(e, t = f) {
		if (super(), Object.defineProperties(this, {
			_intern: { value: /* @__PURE__ */ new Map() },
			_key: { value: t }
		}), e != null) for (let [t, n] of e) this.set(t, n);
	}
	get(e) {
		return super.get(l(this, e));
	}
	has(e) {
		return super.has(l(this, e));
	}
	set(e, t) {
		return super.set(u(this, e), t);
	}
	delete(e) {
		return super.delete(d(this, e));
	}
};
function l({ _intern: e, _key: t }, n) {
	let r = t(n);
	return e.has(r) ? e.get(r) : n;
}
function u({ _intern: e, _key: t }, n) {
	let r = t(n);
	return e.has(r) ? e.get(r) : (e.set(r, n), n);
}
function d({ _intern: e, _key: t }, n) {
	let r = t(n);
	return e.has(r) && (n = e.get(r), e.delete(r)), n;
}
function f(e) {
	return typeof e == "object" && e ? e.valueOf() : e;
}
//#endregion
//#region node_modules/d3-array/src/ticks.js
var p = Math.sqrt(50), m = Math.sqrt(10), h = Math.sqrt(2);
function g(e, t, n) {
	let r = (t - e) / Math.max(0, n), i = Math.floor(Math.log10(r)), a = r / 10 ** i, o = a >= p ? 10 : a >= m ? 5 : a >= h ? 2 : 1, s, c, l;
	return i < 0 ? (l = 10 ** -i / o, s = Math.round(e * l), c = Math.round(t * l), s / l < e && ++s, c / l > t && --c, l = -l) : (l = 10 ** i * o, s = Math.round(e / l), c = Math.round(t / l), s * l < e && ++s, c * l > t && --c), c < s && .5 <= n && n < 2 ? g(e, t, n * 2) : [
		s,
		c,
		l
	];
}
function _(e, t, n) {
	if (t = +t, e = +e, n = +n, !(n > 0)) return [];
	if (e === t) return [e];
	let r = t < e, [i, a, o] = r ? g(t, e, n) : g(e, t, n);
	if (!(a >= i)) return [];
	let s = a - i + 1, c = Array(s);
	if (r) {
		if (o < 0) for (let e = 0; e < s; ++e) c[e] = (a - e) / -o;
		else for (let e = 0; e < s; ++e) c[e] = (a - e) * o;
	} else if (o < 0) for (let e = 0; e < s; ++e) c[e] = (i + e) / -o;
	else for (let e = 0; e < s; ++e) c[e] = (i + e) * o;
	return c;
}
function v(e, t, n) {
	return t = +t, e = +e, n = +n, g(e, t, n)[2];
}
function y(e, t, n) {
	t = +t, e = +e, n = +n;
	let r = t < e, i = r ? v(t, e, n) : v(e, t, n);
	return (r ? -1 : 1) * (i < 0 ? 1 / -i : i);
}
//#endregion
//#region node_modules/d3-array/src/max.js
function b(e, t) {
	let n;
	if (t === void 0) for (let t of e) t != null && (n < t || n === void 0 && t >= t) && (n = t);
	else {
		let r = -1;
		for (let i of e) (i = t(i, ++r, e)) != null && (n < i || n === void 0 && i >= i) && (n = i);
	}
	return n;
}
//#endregion
//#region node_modules/d3-array/src/range.js
function x(e, t, n) {
	e = +e, t = +t, n = (i = arguments.length) < 2 ? (t = e, e = 0, 1) : i < 3 ? 1 : +n;
	for (var r = -1, i = Math.max(0, Math.ceil((t - e) / n)) | 0, a = Array(i); ++r < i;) a[r] = e + r * n;
	return a;
}
//#endregion
//#region node_modules/d3-axis/src/identity.js
function S(e) {
	return e;
}
//#endregion
//#region node_modules/d3-axis/src/axis.js
var C = 1, w = 2, T = 3, E = 4, D = 1e-6;
function O(e) {
	return "translate(" + e + ",0)";
}
function k(e) {
	return "translate(0," + e + ")";
}
function A(e) {
	return (t) => +e(t);
}
function ee(e, t) {
	return t = Math.max(0, e.bandwidth() - t * 2) / 2, e.round() && (t = Math.round(t)), (n) => +e(n) + t;
}
function te() {
	return !this.__axis;
}
function j(e, t) {
	var n = [], r = null, i = null, a = 6, o = 6, s = 3, c = typeof window < "u" && window.devicePixelRatio > 1 ? 0 : .5, l = e === C || e === E ? -1 : 1, u = e === E || e === w ? "x" : "y", d = e === C || e === T ? O : k;
	function f(f) {
		var p = r ?? (t.ticks ? t.ticks.apply(t, n) : t.domain()), m = i ?? (t.tickFormat ? t.tickFormat.apply(t, n) : S), h = Math.max(a, 0) + s, g = t.range(), _ = +g[0] + c, v = +g[g.length - 1] + c, y = (t.bandwidth ? ee : A)(t.copy(), c), b = f.selection ? f.selection() : f, x = b.selectAll(".domain").data([null]), O = b.selectAll(".tick").data(p, t).order(), k = O.exit(), j = O.enter().append("g").attr("class", "tick"), M = O.select("line"), N = O.select("text");
		x = x.merge(x.enter().insert("path", ".tick").attr("class", "domain").attr("stroke", "currentColor")), O = O.merge(j), M = M.merge(j.append("line").attr("stroke", "currentColor").attr(u + "2", l * a)), N = N.merge(j.append("text").attr("fill", "currentColor").attr(u, l * h).attr("dy", e === C ? "0em" : e === T ? "0.71em" : "0.32em")), f !== b && (x = x.transition(f), O = O.transition(f), M = M.transition(f), N = N.transition(f), k = k.transition(f).attr("opacity", D).attr("transform", function(e) {
			return isFinite(e = y(e)) ? d(e + c) : this.getAttribute("transform");
		}), j.attr("opacity", D).attr("transform", function(e) {
			var t = this.parentNode.__axis;
			return d((t && isFinite(t = t(e)) ? t : y(e)) + c);
		})), k.remove(), x.attr("d", e === E || e === w ? o ? "M" + l * o + "," + _ + "H" + c + "V" + v + "H" + l * o : "M" + c + "," + _ + "V" + v : o ? "M" + _ + "," + l * o + "V" + c + "H" + v + "V" + l * o : "M" + _ + "," + c + "H" + v), O.attr("opacity", 1).attr("transform", function(e) {
			return d(y(e) + c);
		}), M.attr(u + "2", l * a), N.attr(u, l * h).text(m), b.filter(te).attr("fill", "none").attr("font-size", 10).attr("font-family", "sans-serif").attr("text-anchor", e === w ? "start" : e === E ? "end" : "middle"), b.each(function() {
			this.__axis = y;
		});
	}
	return f.scale = function(e) {
		return arguments.length ? (t = e, f) : t;
	}, f.ticks = function() {
		return n = Array.from(arguments), f;
	}, f.tickArguments = function(e) {
		return arguments.length ? (n = e == null ? [] : Array.from(e), f) : n.slice();
	}, f.tickValues = function(e) {
		return arguments.length ? (r = e == null ? null : Array.from(e), f) : r && r.slice();
	}, f.tickFormat = function(e) {
		return arguments.length ? (i = e, f) : i;
	}, f.tickSize = function(e) {
		return arguments.length ? (a = o = +e, f) : a;
	}, f.tickSizeInner = function(e) {
		return arguments.length ? (a = +e, f) : a;
	}, f.tickSizeOuter = function(e) {
		return arguments.length ? (o = +e, f) : o;
	}, f.tickPadding = function(e) {
		return arguments.length ? (s = +e, f) : s;
	}, f.offset = function(e) {
		return arguments.length ? (c = +e, f) : c;
	}, f;
}
function M(e) {
	return j(T, e);
}
//#endregion
//#region node_modules/d3-dispatch/src/dispatch.js
var N = { value: () => {} };
function ne() {
	for (var e = 0, t = arguments.length, n = {}, r; e < t; ++e) {
		if (!(r = arguments[e] + "") || r in n || /[\s.]/.test(r)) throw Error("illegal type: " + r);
		n[r] = [];
	}
	return new P(n);
}
function P(e) {
	this._ = e;
}
function re(e, t) {
	return e.trim().split(/^|\s+/).map(function(e) {
		var n = "", r = e.indexOf(".");
		if (r >= 0 && (n = e.slice(r + 1), e = e.slice(0, r)), e && !t.hasOwnProperty(e)) throw Error("unknown type: " + e);
		return {
			type: e,
			name: n
		};
	});
}
P.prototype = ne.prototype = {
	constructor: P,
	on: function(e, t) {
		var n = this._, r = re(e + "", n), i, a = -1, o = r.length;
		if (arguments.length < 2) {
			for (; ++a < o;) if ((i = (e = r[a]).type) && (i = ie(n[i], e.name))) return i;
			return;
		}
		if (t != null && typeof t != "function") throw Error("invalid callback: " + t);
		for (; ++a < o;) if (i = (e = r[a]).type) n[i] = F(n[i], e.name, t);
		else if (t == null) for (i in n) n[i] = F(n[i], e.name, null);
		return this;
	},
	copy: function() {
		var e = {}, t = this._;
		for (var n in t) e[n] = t[n].slice();
		return new P(e);
	},
	call: function(e, t) {
		if ((i = arguments.length - 2) > 0) for (var n = Array(i), r = 0, i, a; r < i; ++r) n[r] = arguments[r + 2];
		if (!this._.hasOwnProperty(e)) throw Error("unknown type: " + e);
		for (a = this._[e], r = 0, i = a.length; r < i; ++r) a[r].value.apply(t, n);
	},
	apply: function(e, t, n) {
		if (!this._.hasOwnProperty(e)) throw Error("unknown type: " + e);
		for (var r = this._[e], i = 0, a = r.length; i < a; ++i) r[i].value.apply(t, n);
	}
};
function ie(e, t) {
	for (var n = 0, r = e.length, i; n < r; ++n) if ((i = e[n]).name === t) return i.value;
}
function F(e, t, n) {
	for (var r = 0, i = e.length; r < i; ++r) if (e[r].name === t) {
		e[r] = N, e = e.slice(0, r).concat(e.slice(r + 1));
		break;
	}
	return n != null && e.push({
		name: t,
		value: n
	}), e;
}
var ae = {
	svg: "http://www.w3.org/2000/svg",
	xhtml: "http://www.w3.org/1999/xhtml",
	xlink: "http://www.w3.org/1999/xlink",
	xml: "http://www.w3.org/XML/1998/namespace",
	xmlns: "http://www.w3.org/2000/xmlns/"
};
//#endregion
//#region node_modules/d3-selection/src/namespace.js
function I(e) {
	var t = e += "", n = t.indexOf(":");
	return n >= 0 && (t = e.slice(0, n)) !== "xmlns" && (e = e.slice(n + 1)), ae.hasOwnProperty(t) ? {
		space: ae[t],
		local: e
	} : e;
}
//#endregion
//#region node_modules/d3-selection/src/creator.js
function oe(e) {
	return function() {
		var t = this.ownerDocument, n = this.namespaceURI;
		return n === "http://www.w3.org/1999/xhtml" && t.documentElement.namespaceURI === "http://www.w3.org/1999/xhtml" ? t.createElement(e) : t.createElementNS(n, e);
	};
}
function se(e) {
	return function() {
		return this.ownerDocument.createElementNS(e.space, e.local);
	};
}
function ce(e) {
	var t = I(e);
	return (t.local ? se : oe)(t);
}
//#endregion
//#region node_modules/d3-selection/src/selector.js
function le() {}
function ue(e) {
	return e == null ? le : function() {
		return this.querySelector(e);
	};
}
//#endregion
//#region node_modules/d3-selection/src/selection/select.js
function de(e) {
	typeof e != "function" && (e = ue(e));
	for (var t = this._groups, n = t.length, r = Array(n), i = 0; i < n; ++i) for (var a = t[i], o = a.length, s = r[i] = Array(o), c, l, u = 0; u < o; ++u) (c = a[u]) && (l = e.call(c, c.__data__, u, a)) && ("__data__" in c && (l.__data__ = c.__data__), s[u] = l);
	return new L(r, this._parents);
}
//#endregion
//#region node_modules/d3-selection/src/array.js
function fe(e) {
	return e == null ? [] : Array.isArray(e) ? e : Array.from(e);
}
//#endregion
//#region node_modules/d3-selection/src/selectorAll.js
function pe() {
	return [];
}
function me(e) {
	return e == null ? pe : function() {
		return this.querySelectorAll(e);
	};
}
//#endregion
//#region node_modules/d3-selection/src/selection/selectAll.js
function he(e) {
	return function() {
		return fe(e.apply(this, arguments));
	};
}
function ge(e) {
	e = typeof e == "function" ? he(e) : me(e);
	for (var t = this._groups, n = t.length, r = [], i = [], a = 0; a < n; ++a) for (var o = t[a], s = o.length, c, l = 0; l < s; ++l) (c = o[l]) && (r.push(e.call(c, c.__data__, l, o)), i.push(c));
	return new L(r, i);
}
//#endregion
//#region node_modules/d3-selection/src/matcher.js
function _e(e) {
	return function() {
		return this.matches(e);
	};
}
function ve(e) {
	return function(t) {
		return t.matches(e);
	};
}
//#endregion
//#region node_modules/d3-selection/src/selection/selectChild.js
var ye = Array.prototype.find;
function be(e) {
	return function() {
		return ye.call(this.children, e);
	};
}
function xe() {
	return this.firstElementChild;
}
function Se(e) {
	return this.select(e == null ? xe : be(typeof e == "function" ? e : ve(e)));
}
//#endregion
//#region node_modules/d3-selection/src/selection/selectChildren.js
var Ce = Array.prototype.filter;
function we() {
	return Array.from(this.children);
}
function Te(e) {
	return function() {
		return Ce.call(this.children, e);
	};
}
function Ee(e) {
	return this.selectAll(e == null ? we : Te(typeof e == "function" ? e : ve(e)));
}
//#endregion
//#region node_modules/d3-selection/src/selection/filter.js
function De(e) {
	typeof e != "function" && (e = _e(e));
	for (var t = this._groups, n = t.length, r = Array(n), i = 0; i < n; ++i) for (var a = t[i], o = a.length, s = r[i] = [], c, l = 0; l < o; ++l) (c = a[l]) && e.call(c, c.__data__, l, a) && s.push(c);
	return new L(r, this._parents);
}
//#endregion
//#region node_modules/d3-selection/src/selection/sparse.js
function Oe(e) {
	return Array(e.length);
}
//#endregion
//#region node_modules/d3-selection/src/selection/enter.js
function ke() {
	return new L(this._enter || this._groups.map(Oe), this._parents);
}
function Ae(e, t) {
	this.ownerDocument = e.ownerDocument, this.namespaceURI = e.namespaceURI, this._next = null, this._parent = e, this.__data__ = t;
}
Ae.prototype = {
	constructor: Ae,
	appendChild: function(e) {
		return this._parent.insertBefore(e, this._next);
	},
	insertBefore: function(e, t) {
		return this._parent.insertBefore(e, t);
	},
	querySelector: function(e) {
		return this._parent.querySelector(e);
	},
	querySelectorAll: function(e) {
		return this._parent.querySelectorAll(e);
	}
};
//#endregion
//#region node_modules/d3-selection/src/constant.js
function je(e) {
	return function() {
		return e;
	};
}
//#endregion
//#region node_modules/d3-selection/src/selection/data.js
function Me(e, t, n, r, i, a) {
	for (var o = 0, s, c = t.length, l = a.length; o < l; ++o) (s = t[o]) ? (s.__data__ = a[o], r[o] = s) : n[o] = new Ae(e, a[o]);
	for (; o < c; ++o) (s = t[o]) && (i[o] = s);
}
function Ne(e, t, n, r, i, a, o) {
	var s, c, l = /* @__PURE__ */ new Map(), u = t.length, d = a.length, f = Array(u), p;
	for (s = 0; s < u; ++s) (c = t[s]) && (f[s] = p = o.call(c, c.__data__, s, t) + "", l.has(p) ? i[s] = c : l.set(p, c));
	for (s = 0; s < d; ++s) p = o.call(e, a[s], s, a) + "", (c = l.get(p)) ? (r[s] = c, c.__data__ = a[s], l.delete(p)) : n[s] = new Ae(e, a[s]);
	for (s = 0; s < u; ++s) (c = t[s]) && l.get(f[s]) === c && (i[s] = c);
}
function Pe(e) {
	return e.__data__;
}
function Fe(e, t) {
	if (!arguments.length) return Array.from(this, Pe);
	var n = t ? Ne : Me, r = this._parents, i = this._groups;
	typeof e != "function" && (e = je(e));
	for (var a = i.length, o = Array(a), s = Array(a), c = Array(a), l = 0; l < a; ++l) {
		var u = r[l], d = i[l], f = d.length, p = Ie(e.call(u, u && u.__data__, l, r)), m = p.length, h = s[l] = Array(m), g = o[l] = Array(m);
		n(u, d, h, g, c[l] = Array(f), p, t);
		for (var _ = 0, v = 0, y, b; _ < m; ++_) if (y = h[_]) {
			for (_ >= v && (v = _ + 1); !(b = g[v]) && ++v < m;);
			y._next = b || null;
		}
	}
	return o = new L(o, r), o._enter = s, o._exit = c, o;
}
function Ie(e) {
	return typeof e == "object" && "length" in e ? e : Array.from(e);
}
//#endregion
//#region node_modules/d3-selection/src/selection/exit.js
function Le() {
	return new L(this._exit || this._groups.map(Oe), this._parents);
}
//#endregion
//#region node_modules/d3-selection/src/selection/join.js
function Re(e, t, n) {
	var r = this.enter(), i = this, a = this.exit();
	return typeof e == "function" ? (r = e(r), r &&= r.selection()) : r = r.append(e + ""), t != null && (i = t(i), i &&= i.selection()), n == null ? a.remove() : n(a), r && i ? r.merge(i).order() : i;
}
//#endregion
//#region node_modules/d3-selection/src/selection/merge.js
function ze(e) {
	for (var t = e.selection ? e.selection() : e, n = this._groups, r = t._groups, i = n.length, a = r.length, o = Math.min(i, a), s = Array(i), c = 0; c < o; ++c) for (var l = n[c], u = r[c], d = l.length, f = s[c] = Array(d), p, m = 0; m < d; ++m) (p = l[m] || u[m]) && (f[m] = p);
	for (; c < i; ++c) s[c] = n[c];
	return new L(s, this._parents);
}
//#endregion
//#region node_modules/d3-selection/src/selection/order.js
function Be() {
	for (var e = this._groups, t = -1, n = e.length; ++t < n;) for (var r = e[t], i = r.length - 1, a = r[i], o; --i >= 0;) (o = r[i]) && (a && o.compareDocumentPosition(a) ^ 4 && a.parentNode.insertBefore(o, a), a = o);
	return this;
}
//#endregion
//#region node_modules/d3-selection/src/selection/sort.js
function Ve(e) {
	e ||= He;
	function t(t, n) {
		return t && n ? e(t.__data__, n.__data__) : !t - !n;
	}
	for (var n = this._groups, r = n.length, i = Array(r), a = 0; a < r; ++a) {
		for (var o = n[a], s = o.length, c = i[a] = Array(s), l, u = 0; u < s; ++u) (l = o[u]) && (c[u] = l);
		c.sort(t);
	}
	return new L(i, this._parents).order();
}
function He(e, t) {
	return e < t ? -1 : e > t ? 1 : e >= t ? 0 : NaN;
}
//#endregion
//#region node_modules/d3-selection/src/selection/call.js
function Ue() {
	var e = arguments[0];
	return arguments[0] = this, e.apply(null, arguments), this;
}
//#endregion
//#region node_modules/d3-selection/src/selection/nodes.js
function We() {
	return Array.from(this);
}
//#endregion
//#region node_modules/d3-selection/src/selection/node.js
function Ge() {
	for (var e = this._groups, t = 0, n = e.length; t < n; ++t) for (var r = e[t], i = 0, a = r.length; i < a; ++i) {
		var o = r[i];
		if (o) return o;
	}
	return null;
}
//#endregion
//#region node_modules/d3-selection/src/selection/size.js
function Ke() {
	let e = 0;
	for (let t of this) ++e;
	return e;
}
//#endregion
//#region node_modules/d3-selection/src/selection/empty.js
function qe() {
	return !this.node();
}
//#endregion
//#region node_modules/d3-selection/src/selection/each.js
function Je(e) {
	for (var t = this._groups, n = 0, r = t.length; n < r; ++n) for (var i = t[n], a = 0, o = i.length, s; a < o; ++a) (s = i[a]) && e.call(s, s.__data__, a, i);
	return this;
}
//#endregion
//#region node_modules/d3-selection/src/selection/attr.js
function Ye(e) {
	return function() {
		this.removeAttribute(e);
	};
}
function Xe(e) {
	return function() {
		this.removeAttributeNS(e.space, e.local);
	};
}
function Ze(e, t) {
	return function() {
		this.setAttribute(e, t);
	};
}
function Qe(e, t) {
	return function() {
		this.setAttributeNS(e.space, e.local, t);
	};
}
function $e(e, t) {
	return function() {
		var n = t.apply(this, arguments);
		n == null ? this.removeAttribute(e) : this.setAttribute(e, n);
	};
}
function et(e, t) {
	return function() {
		var n = t.apply(this, arguments);
		n == null ? this.removeAttributeNS(e.space, e.local) : this.setAttributeNS(e.space, e.local, n);
	};
}
function tt(e, t) {
	var n = I(e);
	if (arguments.length < 2) {
		var r = this.node();
		return n.local ? r.getAttributeNS(n.space, n.local) : r.getAttribute(n);
	}
	return this.each((t == null ? n.local ? Xe : Ye : typeof t == "function" ? n.local ? et : $e : n.local ? Qe : Ze)(n, t));
}
//#endregion
//#region node_modules/d3-selection/src/window.js
function nt(e) {
	return e.ownerDocument && e.ownerDocument.defaultView || e.document && e || e.defaultView;
}
//#endregion
//#region node_modules/d3-selection/src/selection/style.js
function rt(e) {
	return function() {
		this.style.removeProperty(e);
	};
}
function it(e, t, n) {
	return function() {
		this.style.setProperty(e, t, n);
	};
}
function at(e, t, n) {
	return function() {
		var r = t.apply(this, arguments);
		r == null ? this.style.removeProperty(e) : this.style.setProperty(e, r, n);
	};
}
function ot(e, t, n) {
	return arguments.length > 1 ? this.each((t == null ? rt : typeof t == "function" ? at : it)(e, t, n ?? "")) : st(this.node(), e);
}
function st(e, t) {
	return e.style.getPropertyValue(t) || nt(e).getComputedStyle(e, null).getPropertyValue(t);
}
//#endregion
//#region node_modules/d3-selection/src/selection/property.js
function ct(e) {
	return function() {
		delete this[e];
	};
}
function lt(e, t) {
	return function() {
		this[e] = t;
	};
}
function ut(e, t) {
	return function() {
		var n = t.apply(this, arguments);
		n == null ? delete this[e] : this[e] = n;
	};
}
function dt(e, t) {
	return arguments.length > 1 ? this.each((t == null ? ct : typeof t == "function" ? ut : lt)(e, t)) : this.node()[e];
}
//#endregion
//#region node_modules/d3-selection/src/selection/classed.js
function ft(e) {
	return e.trim().split(/^|\s+/);
}
function pt(e) {
	return e.classList || new mt(e);
}
function mt(e) {
	this._node = e, this._names = ft(e.getAttribute("class") || "");
}
mt.prototype = {
	add: function(e) {
		this._names.indexOf(e) < 0 && (this._names.push(e), this._node.setAttribute("class", this._names.join(" ")));
	},
	remove: function(e) {
		var t = this._names.indexOf(e);
		t >= 0 && (this._names.splice(t, 1), this._node.setAttribute("class", this._names.join(" ")));
	},
	contains: function(e) {
		return this._names.indexOf(e) >= 0;
	}
};
function ht(e, t) {
	for (var n = pt(e), r = -1, i = t.length; ++r < i;) n.add(t[r]);
}
function gt(e, t) {
	for (var n = pt(e), r = -1, i = t.length; ++r < i;) n.remove(t[r]);
}
function _t(e) {
	return function() {
		ht(this, e);
	};
}
function vt(e) {
	return function() {
		gt(this, e);
	};
}
function yt(e, t) {
	return function() {
		(t.apply(this, arguments) ? ht : gt)(this, e);
	};
}
function bt(e, t) {
	var n = ft(e + "");
	if (arguments.length < 2) {
		for (var r = pt(this.node()), i = -1, a = n.length; ++i < a;) if (!r.contains(n[i])) return !1;
		return !0;
	}
	return this.each((typeof t == "function" ? yt : t ? _t : vt)(n, t));
}
//#endregion
//#region node_modules/d3-selection/src/selection/text.js
function xt() {
	this.textContent = "";
}
function St(e) {
	return function() {
		this.textContent = e;
	};
}
function Ct(e) {
	return function() {
		var t = e.apply(this, arguments);
		this.textContent = t ?? "";
	};
}
function wt(e) {
	return arguments.length ? this.each(e == null ? xt : (typeof e == "function" ? Ct : St)(e)) : this.node().textContent;
}
//#endregion
//#region node_modules/d3-selection/src/selection/html.js
function Tt() {
	this.innerHTML = "";
}
function Et(e) {
	return function() {
		this.innerHTML = e;
	};
}
function Dt(e) {
	return function() {
		var t = e.apply(this, arguments);
		this.innerHTML = t ?? "";
	};
}
function Ot(e) {
	return arguments.length ? this.each(e == null ? Tt : (typeof e == "function" ? Dt : Et)(e)) : this.node().innerHTML;
}
//#endregion
//#region node_modules/d3-selection/src/selection/raise.js
function kt() {
	this.nextSibling && this.parentNode.appendChild(this);
}
function At() {
	return this.each(kt);
}
//#endregion
//#region node_modules/d3-selection/src/selection/lower.js
function jt() {
	this.previousSibling && this.parentNode.insertBefore(this, this.parentNode.firstChild);
}
function Mt() {
	return this.each(jt);
}
//#endregion
//#region node_modules/d3-selection/src/selection/append.js
function Nt(e) {
	var t = typeof e == "function" ? e : ce(e);
	return this.select(function() {
		return this.appendChild(t.apply(this, arguments));
	});
}
//#endregion
//#region node_modules/d3-selection/src/selection/insert.js
function Pt() {
	return null;
}
function Ft(e, t) {
	var n = typeof e == "function" ? e : ce(e), r = t == null ? Pt : typeof t == "function" ? t : ue(t);
	return this.select(function() {
		return this.insertBefore(n.apply(this, arguments), r.apply(this, arguments) || null);
	});
}
//#endregion
//#region node_modules/d3-selection/src/selection/remove.js
function It() {
	var e = this.parentNode;
	e && e.removeChild(this);
}
function Lt() {
	return this.each(It);
}
//#endregion
//#region node_modules/d3-selection/src/selection/clone.js
function Rt() {
	var e = this.cloneNode(!1), t = this.parentNode;
	return t ? t.insertBefore(e, this.nextSibling) : e;
}
function zt() {
	var e = this.cloneNode(!0), t = this.parentNode;
	return t ? t.insertBefore(e, this.nextSibling) : e;
}
function Bt(e) {
	return this.select(e ? zt : Rt);
}
//#endregion
//#region node_modules/d3-selection/src/selection/datum.js
function Vt(e) {
	return arguments.length ? this.property("__data__", e) : this.node().__data__;
}
//#endregion
//#region node_modules/d3-selection/src/selection/on.js
function Ht(e) {
	return function(t) {
		e.call(this, t, this.__data__);
	};
}
function Ut(e) {
	return e.trim().split(/^|\s+/).map(function(e) {
		var t = "", n = e.indexOf(".");
		return n >= 0 && (t = e.slice(n + 1), e = e.slice(0, n)), {
			type: e,
			name: t
		};
	});
}
function Wt(e) {
	return function() {
		var t = this.__on;
		if (t) {
			for (var n = 0, r = -1, i = t.length, a; n < i; ++n) a = t[n], (!e.type || a.type === e.type) && a.name === e.name ? this.removeEventListener(a.type, a.listener, a.options) : t[++r] = a;
			++r ? t.length = r : delete this.__on;
		}
	};
}
function Gt(e, t, n) {
	return function() {
		var r = this.__on, i, a = Ht(t);
		if (r) {
			for (var o = 0, s = r.length; o < s; ++o) if ((i = r[o]).type === e.type && i.name === e.name) {
				this.removeEventListener(i.type, i.listener, i.options), this.addEventListener(i.type, i.listener = a, i.options = n), i.value = t;
				return;
			}
		}
		this.addEventListener(e.type, a, n), i = {
			type: e.type,
			name: e.name,
			value: t,
			listener: a,
			options: n
		}, r ? r.push(i) : this.__on = [i];
	};
}
function Kt(e, t, n) {
	var r = Ut(e + ""), i, a = r.length, o;
	if (arguments.length < 2) {
		var s = this.node().__on;
		if (s) {
			for (var c = 0, l = s.length, u; c < l; ++c) for (i = 0, u = s[c]; i < a; ++i) if ((o = r[i]).type === u.type && o.name === u.name) return u.value;
		}
		return;
	}
	for (s = t ? Gt : Wt, i = 0; i < a; ++i) this.each(s(r[i], t, n));
	return this;
}
//#endregion
//#region node_modules/d3-selection/src/selection/dispatch.js
function qt(e, t, n) {
	var r = nt(e), i = r.CustomEvent;
	typeof i == "function" ? i = new i(t, n) : (i = r.document.createEvent("Event"), n ? (i.initEvent(t, n.bubbles, n.cancelable), i.detail = n.detail) : i.initEvent(t, !1, !1)), e.dispatchEvent(i);
}
function Jt(e, t) {
	return function() {
		return qt(this, e, t);
	};
}
function Yt(e, t) {
	return function() {
		return qt(this, e, t.apply(this, arguments));
	};
}
function Xt(e, t) {
	return this.each((typeof t == "function" ? Yt : Jt)(e, t));
}
//#endregion
//#region node_modules/d3-selection/src/selection/iterator.js
function* Zt() {
	for (var e = this._groups, t = 0, n = e.length; t < n; ++t) for (var r = e[t], i = 0, a = r.length, o; i < a; ++i) (o = r[i]) && (yield o);
}
//#endregion
//#region node_modules/d3-selection/src/selection/index.js
var Qt = [null];
function L(e, t) {
	this._groups = e, this._parents = t;
}
function $t() {
	return new L([[document.documentElement]], Qt);
}
function en() {
	return this;
}
L.prototype = $t.prototype = {
	constructor: L,
	select: de,
	selectAll: ge,
	selectChild: Se,
	selectChildren: Ee,
	filter: De,
	data: Fe,
	enter: ke,
	exit: Le,
	join: Re,
	merge: ze,
	selection: en,
	order: Be,
	sort: Ve,
	call: Ue,
	nodes: We,
	node: Ge,
	size: Ke,
	empty: qe,
	each: Je,
	attr: tt,
	style: ot,
	property: dt,
	classed: bt,
	text: wt,
	html: Ot,
	raise: At,
	lower: Mt,
	append: Nt,
	insert: Ft,
	remove: Lt,
	clone: Bt,
	datum: Vt,
	on: Kt,
	dispatch: Xt,
	[Symbol.iterator]: Zt
};
//#endregion
//#region node_modules/d3-selection/src/select.js
function R(e) {
	return typeof e == "string" ? new L([[document.querySelector(e)]], [document.documentElement]) : new L([[e]], Qt);
}
//#endregion
//#region node_modules/d3-selection/src/sourceEvent.js
function tn(e) {
	let t;
	for (; t = e.sourceEvent;) e = t;
	return e;
}
//#endregion
//#region node_modules/d3-selection/src/pointer.js
function nn(e, t) {
	if (e = tn(e), t === void 0 && (t = e.currentTarget), t) {
		var n = t.ownerSVGElement || t;
		if (n.createSVGPoint) {
			var r = n.createSVGPoint();
			return r.x = e.clientX, r.y = e.clientY, r = r.matrixTransform(t.getScreenCTM().inverse()), [r.x, r.y];
		}
		if (t.getBoundingClientRect) {
			var i = t.getBoundingClientRect();
			return [e.clientX - i.left - t.clientLeft, e.clientY - i.top - t.clientTop];
		}
	}
	return [e.pageX, e.pageY];
}
//#endregion
//#region node_modules/d3-drag/src/noevent.js
var rn = {
	capture: !0,
	passive: !1
};
function an(e) {
	e.preventDefault(), e.stopImmediatePropagation();
}
//#endregion
//#region node_modules/d3-drag/src/nodrag.js
function on(e) {
	var t = e.document.documentElement, n = R(e).on("dragstart.drag", an, rn);
	"onselectstart" in t ? n.on("selectstart.drag", an, rn) : (t.__noselect = t.style.MozUserSelect, t.style.MozUserSelect = "none");
}
function sn(e, t) {
	var n = e.document.documentElement, r = R(e).on("dragstart.drag", null);
	t && (r.on("click.drag", an, rn), setTimeout(function() {
		r.on("click.drag", null);
	}, 0)), "onselectstart" in n ? r.on("selectstart.drag", null) : (n.style.MozUserSelect = n.__noselect, delete n.__noselect);
}
//#endregion
//#region node_modules/d3-color/src/define.js
function cn(e, t, n) {
	e.prototype = t.prototype = n, n.constructor = e;
}
function ln(e, t) {
	var n = Object.create(e.prototype);
	for (var r in t) n[r] = t[r];
	return n;
}
//#endregion
//#region node_modules/d3-color/src/color.js
function un() {}
var dn = .7, fn = 1 / dn, pn = "\\s*([+-]?\\d+)\\s*", mn = "\\s*([+-]?(?:\\d*\\.)?\\d+(?:[eE][+-]?\\d+)?)\\s*", z = "\\s*([+-]?(?:\\d*\\.)?\\d+(?:[eE][+-]?\\d+)?)%\\s*", hn = /^#([0-9a-f]{3,8})$/, gn = RegExp(`^rgb\\(${pn},${pn},${pn}\\)$`), _n = RegExp(`^rgb\\(${z},${z},${z}\\)$`), vn = RegExp(`^rgba\\(${pn},${pn},${pn},${mn}\\)$`), yn = RegExp(`^rgba\\(${z},${z},${z},${mn}\\)$`), bn = RegExp(`^hsl\\(${mn},${z},${z}\\)$`), xn = RegExp(`^hsla\\(${mn},${z},${z},${mn}\\)$`), Sn = {
	aliceblue: 15792383,
	antiquewhite: 16444375,
	aqua: 65535,
	aquamarine: 8388564,
	azure: 15794175,
	beige: 16119260,
	bisque: 16770244,
	black: 0,
	blanchedalmond: 16772045,
	blue: 255,
	blueviolet: 9055202,
	brown: 10824234,
	burlywood: 14596231,
	cadetblue: 6266528,
	chartreuse: 8388352,
	chocolate: 13789470,
	coral: 16744272,
	cornflowerblue: 6591981,
	cornsilk: 16775388,
	crimson: 14423100,
	cyan: 65535,
	darkblue: 139,
	darkcyan: 35723,
	darkgoldenrod: 12092939,
	darkgray: 11119017,
	darkgreen: 25600,
	darkgrey: 11119017,
	darkkhaki: 12433259,
	darkmagenta: 9109643,
	darkolivegreen: 5597999,
	darkorange: 16747520,
	darkorchid: 10040012,
	darkred: 9109504,
	darksalmon: 15308410,
	darkseagreen: 9419919,
	darkslateblue: 4734347,
	darkslategray: 3100495,
	darkslategrey: 3100495,
	darkturquoise: 52945,
	darkviolet: 9699539,
	deeppink: 16716947,
	deepskyblue: 49151,
	dimgray: 6908265,
	dimgrey: 6908265,
	dodgerblue: 2003199,
	firebrick: 11674146,
	floralwhite: 16775920,
	forestgreen: 2263842,
	fuchsia: 16711935,
	gainsboro: 14474460,
	ghostwhite: 16316671,
	gold: 16766720,
	goldenrod: 14329120,
	gray: 8421504,
	green: 32768,
	greenyellow: 11403055,
	grey: 8421504,
	honeydew: 15794160,
	hotpink: 16738740,
	indianred: 13458524,
	indigo: 4915330,
	ivory: 16777200,
	khaki: 15787660,
	lavender: 15132410,
	lavenderblush: 16773365,
	lawngreen: 8190976,
	lemonchiffon: 16775885,
	lightblue: 11393254,
	lightcoral: 15761536,
	lightcyan: 14745599,
	lightgoldenrodyellow: 16448210,
	lightgray: 13882323,
	lightgreen: 9498256,
	lightgrey: 13882323,
	lightpink: 16758465,
	lightsalmon: 16752762,
	lightseagreen: 2142890,
	lightskyblue: 8900346,
	lightslategray: 7833753,
	lightslategrey: 7833753,
	lightsteelblue: 11584734,
	lightyellow: 16777184,
	lime: 65280,
	limegreen: 3329330,
	linen: 16445670,
	magenta: 16711935,
	maroon: 8388608,
	mediumaquamarine: 6737322,
	mediumblue: 205,
	mediumorchid: 12211667,
	mediumpurple: 9662683,
	mediumseagreen: 3978097,
	mediumslateblue: 8087790,
	mediumspringgreen: 64154,
	mediumturquoise: 4772300,
	mediumvioletred: 13047173,
	midnightblue: 1644912,
	mintcream: 16121850,
	mistyrose: 16770273,
	moccasin: 16770229,
	navajowhite: 16768685,
	navy: 128,
	oldlace: 16643558,
	olive: 8421376,
	olivedrab: 7048739,
	orange: 16753920,
	orangered: 16729344,
	orchid: 14315734,
	palegoldenrod: 15657130,
	palegreen: 10025880,
	paleturquoise: 11529966,
	palevioletred: 14381203,
	papayawhip: 16773077,
	peachpuff: 16767673,
	peru: 13468991,
	pink: 16761035,
	plum: 14524637,
	powderblue: 11591910,
	purple: 8388736,
	rebeccapurple: 6697881,
	red: 16711680,
	rosybrown: 12357519,
	royalblue: 4286945,
	saddlebrown: 9127187,
	salmon: 16416882,
	sandybrown: 16032864,
	seagreen: 3050327,
	seashell: 16774638,
	sienna: 10506797,
	silver: 12632256,
	skyblue: 8900331,
	slateblue: 6970061,
	slategray: 7372944,
	slategrey: 7372944,
	snow: 16775930,
	springgreen: 65407,
	steelblue: 4620980,
	tan: 13808780,
	teal: 32896,
	thistle: 14204888,
	tomato: 16737095,
	turquoise: 4251856,
	violet: 15631086,
	wheat: 16113331,
	white: 16777215,
	whitesmoke: 16119285,
	yellow: 16776960,
	yellowgreen: 10145074
};
cn(un, Dn, {
	copy(e) {
		return Object.assign(new this.constructor(), this, e);
	},
	displayable() {
		return this.rgb().displayable();
	},
	hex: Cn,
	formatHex: Cn,
	formatHex8: wn,
	formatHsl: Tn,
	formatRgb: En,
	toString: En
});
function Cn() {
	return this.rgb().formatHex();
}
function wn() {
	return this.rgb().formatHex8();
}
function Tn() {
	return zn(this).formatHsl();
}
function En() {
	return this.rgb().formatRgb();
}
function Dn(e) {
	var t, n;
	return e = (e + "").trim().toLowerCase(), (t = hn.exec(e)) ? (n = t[1].length, t = parseInt(t[1], 16), n === 6 ? On(t) : n === 3 ? new B(t >> 8 & 15 | t >> 4 & 240, t >> 4 & 15 | t & 240, (t & 15) << 4 | t & 15, 1) : n === 8 ? kn(t >> 24 & 255, t >> 16 & 255, t >> 8 & 255, (t & 255) / 255) : n === 4 ? kn(t >> 12 & 15 | t >> 8 & 240, t >> 8 & 15 | t >> 4 & 240, t >> 4 & 15 | t & 240, ((t & 15) << 4 | t & 15) / 255) : null) : (t = gn.exec(e)) ? new B(t[1], t[2], t[3], 1) : (t = _n.exec(e)) ? new B(t[1] * 255 / 100, t[2] * 255 / 100, t[3] * 255 / 100, 1) : (t = vn.exec(e)) ? kn(t[1], t[2], t[3], t[4]) : (t = yn.exec(e)) ? kn(t[1] * 255 / 100, t[2] * 255 / 100, t[3] * 255 / 100, t[4]) : (t = bn.exec(e)) ? Rn(t[1], t[2] / 100, t[3] / 100, 1) : (t = xn.exec(e)) ? Rn(t[1], t[2] / 100, t[3] / 100, t[4]) : Sn.hasOwnProperty(e) ? On(Sn[e]) : e === "transparent" ? new B(NaN, NaN, NaN, 0) : null;
}
function On(e) {
	return new B(e >> 16 & 255, e >> 8 & 255, e & 255, 1);
}
function kn(e, t, n, r) {
	return r <= 0 && (e = t = n = NaN), new B(e, t, n, r);
}
function An(e) {
	return e instanceof un || (e = Dn(e)), e ? (e = e.rgb(), new B(e.r, e.g, e.b, e.opacity)) : new B();
}
function jn(e, t, n, r) {
	return arguments.length === 1 ? An(e) : new B(e, t, n, r ?? 1);
}
function B(e, t, n, r) {
	this.r = +e, this.g = +t, this.b = +n, this.opacity = +r;
}
cn(B, jn, ln(un, {
	brighter(e) {
		return e = e == null ? fn : fn ** +e, new B(this.r * e, this.g * e, this.b * e, this.opacity);
	},
	darker(e) {
		return e = e == null ? dn : dn ** +e, new B(this.r * e, this.g * e, this.b * e, this.opacity);
	},
	rgb() {
		return this;
	},
	clamp() {
		return new B(In(this.r), In(this.g), In(this.b), Fn(this.opacity));
	},
	displayable() {
		return -.5 <= this.r && this.r < 255.5 && -.5 <= this.g && this.g < 255.5 && -.5 <= this.b && this.b < 255.5 && 0 <= this.opacity && this.opacity <= 1;
	},
	hex: Mn,
	formatHex: Mn,
	formatHex8: Nn,
	formatRgb: Pn,
	toString: Pn
}));
function Mn() {
	return `#${Ln(this.r)}${Ln(this.g)}${Ln(this.b)}`;
}
function Nn() {
	return `#${Ln(this.r)}${Ln(this.g)}${Ln(this.b)}${Ln((isNaN(this.opacity) ? 1 : this.opacity) * 255)}`;
}
function Pn() {
	let e = Fn(this.opacity);
	return `${e === 1 ? "rgb(" : "rgba("}${In(this.r)}, ${In(this.g)}, ${In(this.b)}${e === 1 ? ")" : `, ${e})`}`;
}
function Fn(e) {
	return isNaN(e) ? 1 : Math.max(0, Math.min(1, e));
}
function In(e) {
	return Math.max(0, Math.min(255, Math.round(e) || 0));
}
function Ln(e) {
	return e = In(e), (e < 16 ? "0" : "") + e.toString(16);
}
function Rn(e, t, n, r) {
	return r <= 0 ? e = t = n = NaN : n <= 0 || n >= 1 ? e = t = NaN : t <= 0 && (e = NaN), new V(e, t, n, r);
}
function zn(e) {
	if (e instanceof V) return new V(e.h, e.s, e.l, e.opacity);
	if (e instanceof un || (e = Dn(e)), !e) return new V();
	if (e instanceof V) return e;
	e = e.rgb();
	var t = e.r / 255, n = e.g / 255, r = e.b / 255, i = Math.min(t, n, r), a = Math.max(t, n, r), o = NaN, s = a - i, c = (a + i) / 2;
	return s ? (o = t === a ? (n - r) / s + (n < r) * 6 : n === a ? (r - t) / s + 2 : (t - n) / s + 4, s /= c < .5 ? a + i : 2 - a - i, o *= 60) : s = c > 0 && c < 1 ? 0 : o, new V(o, s, c, e.opacity);
}
function Bn(e, t, n, r) {
	return arguments.length === 1 ? zn(e) : new V(e, t, n, r ?? 1);
}
function V(e, t, n, r) {
	this.h = +e, this.s = +t, this.l = +n, this.opacity = +r;
}
cn(V, Bn, ln(un, {
	brighter(e) {
		return e = e == null ? fn : fn ** +e, new V(this.h, this.s, this.l * e, this.opacity);
	},
	darker(e) {
		return e = e == null ? dn : dn ** +e, new V(this.h, this.s, this.l * e, this.opacity);
	},
	rgb() {
		var e = this.h % 360 + (this.h < 0) * 360, t = isNaN(e) || isNaN(this.s) ? 0 : this.s, n = this.l, r = n + (n < .5 ? n : 1 - n) * t, i = 2 * n - r;
		return new B(Un(e >= 240 ? e - 240 : e + 120, i, r), Un(e, i, r), Un(e < 120 ? e + 240 : e - 120, i, r), this.opacity);
	},
	clamp() {
		return new V(Vn(this.h), Hn(this.s), Hn(this.l), Fn(this.opacity));
	},
	displayable() {
		return (0 <= this.s && this.s <= 1 || isNaN(this.s)) && 0 <= this.l && this.l <= 1 && 0 <= this.opacity && this.opacity <= 1;
	},
	formatHsl() {
		let e = Fn(this.opacity);
		return `${e === 1 ? "hsl(" : "hsla("}${Vn(this.h)}, ${Hn(this.s) * 100}%, ${Hn(this.l) * 100}%${e === 1 ? ")" : `, ${e})`}`;
	}
}));
function Vn(e) {
	return e = (e || 0) % 360, e < 0 ? e + 360 : e;
}
function Hn(e) {
	return Math.max(0, Math.min(1, e || 0));
}
function Un(e, t, n) {
	return (e < 60 ? t + (n - t) * e / 60 : e < 180 ? n : e < 240 ? t + (n - t) * (240 - e) / 60 : t) * 255;
}
//#endregion
//#region node_modules/d3-color/src/math.js
var Wn = Math.PI / 180, Gn = 180 / Math.PI, Kn = 18, qn = .96422, Jn = 1, Yn = .82521, Xn = 4 / 29, Zn = 6 / 29, Qn = 3 * Zn * Zn, $n = Zn * Zn * Zn;
function er(e) {
	if (e instanceof H) return new H(e.l, e.a, e.b, e.opacity);
	if (e instanceof U) return cr(e);
	e instanceof B || (e = An(e));
	var t = ar(e.r), n = ar(e.g), r = ar(e.b), i = nr((.2225045 * t + .7168786 * n + .0606169 * r) / Jn), a, o;
	return t === n && n === r ? a = o = i : (a = nr((.4360747 * t + .3850649 * n + .1430804 * r) / qn), o = nr((.0139322 * t + .0971045 * n + .7141733 * r) / Yn)), new H(116 * i - 16, 500 * (a - i), 200 * (i - o), e.opacity);
}
function tr(e, t, n, r) {
	return arguments.length === 1 ? er(e) : new H(e, t, n, r ?? 1);
}
function H(e, t, n, r) {
	this.l = +e, this.a = +t, this.b = +n, this.opacity = +r;
}
cn(H, tr, ln(un, {
	brighter(e) {
		return new H(this.l + Kn * (e ?? 1), this.a, this.b, this.opacity);
	},
	darker(e) {
		return new H(this.l - Kn * (e ?? 1), this.a, this.b, this.opacity);
	},
	rgb() {
		var e = (this.l + 16) / 116, t = isNaN(this.a) ? e : e + this.a / 500, n = isNaN(this.b) ? e : e - this.b / 200;
		return t = qn * rr(t), e = Jn * rr(e), n = Yn * rr(n), new B(ir(3.1338561 * t - 1.6168667 * e - .4906146 * n), ir(-.9787684 * t + 1.9161415 * e + .033454 * n), ir(.0719453 * t - .2289914 * e + 1.4052427 * n), this.opacity);
	}
}));
function nr(e) {
	return e > $n ? e ** (1 / 3) : e / Qn + Xn;
}
function rr(e) {
	return e > Zn ? e * e * e : Qn * (e - Xn);
}
function ir(e) {
	return 255 * (e <= .0031308 ? 12.92 * e : 1.055 * e ** (1 / 2.4) - .055);
}
function ar(e) {
	return (e /= 255) <= .04045 ? e / 12.92 : ((e + .055) / 1.055) ** 2.4;
}
function or(e) {
	if (e instanceof U) return new U(e.h, e.c, e.l, e.opacity);
	if (e instanceof H || (e = er(e)), e.a === 0 && e.b === 0) return new U(NaN, 0 < e.l && e.l < 100 ? 0 : NaN, e.l, e.opacity);
	var t = Math.atan2(e.b, e.a) * Gn;
	return new U(t < 0 ? t + 360 : t, Math.sqrt(e.a * e.a + e.b * e.b), e.l, e.opacity);
}
function sr(e, t, n, r) {
	return arguments.length === 1 ? or(e) : new U(e, t, n, r ?? 1);
}
function U(e, t, n, r) {
	this.h = +e, this.c = +t, this.l = +n, this.opacity = +r;
}
function cr(e) {
	if (isNaN(e.h)) return new H(e.l, 0, 0, e.opacity);
	var t = e.h * Wn;
	return new H(e.l, Math.cos(t) * e.c, Math.sin(t) * e.c, e.opacity);
}
cn(U, sr, ln(un, {
	brighter(e) {
		return new U(this.h, this.c, this.l + Kn * (e ?? 1), this.opacity);
	},
	darker(e) {
		return new U(this.h, this.c, this.l - Kn * (e ?? 1), this.opacity);
	},
	rgb() {
		return cr(this).rgb();
	}
}));
//#endregion
//#region node_modules/d3-interpolate/src/constant.js
var lr = (e) => () => e;
//#endregion
//#region node_modules/d3-interpolate/src/color.js
function ur(e, t) {
	return function(n) {
		return e + n * t;
	};
}
function dr(e, t, n) {
	return e **= +n, t = t ** +n - e, n = 1 / n, function(r) {
		return (e + r * t) ** +n;
	};
}
function fr(e) {
	return (e = +e) == 1 ? pr : function(t, n) {
		return n - t ? dr(t, n, e) : lr(isNaN(t) ? n : t);
	};
}
function pr(e, t) {
	var n = t - e;
	return n ? ur(e, n) : lr(isNaN(e) ? t : e);
}
//#endregion
//#region node_modules/d3-interpolate/src/rgb.js
var mr = (function e(t) {
	var n = fr(t);
	function r(e, t) {
		var r = n((e = jn(e)).r, (t = jn(t)).r), i = n(e.g, t.g), a = n(e.b, t.b), o = pr(e.opacity, t.opacity);
		return function(t) {
			return e.r = r(t), e.g = i(t), e.b = a(t), e.opacity = o(t), e + "";
		};
	}
	return r.gamma = e, r;
})(1);
//#endregion
//#region node_modules/d3-interpolate/src/numberArray.js
function hr(e, t) {
	t ||= [];
	var n = e ? Math.min(t.length, e.length) : 0, r = t.slice(), i;
	return function(a) {
		for (i = 0; i < n; ++i) r[i] = e[i] * (1 - a) + t[i] * a;
		return r;
	};
}
function gr(e) {
	return ArrayBuffer.isView(e) && !(e instanceof DataView);
}
//#endregion
//#region node_modules/d3-interpolate/src/array.js
function _r(e, t) {
	var n = t ? t.length : 0, r = e ? Math.min(n, e.length) : 0, i = Array(r), a = Array(n), o;
	for (o = 0; o < r; ++o) i[o] = Tr(e[o], t[o]);
	for (; o < n; ++o) a[o] = t[o];
	return function(e) {
		for (o = 0; o < r; ++o) a[o] = i[o](e);
		return a;
	};
}
//#endregion
//#region node_modules/d3-interpolate/src/date.js
function vr(e, t) {
	var n = /* @__PURE__ */ new Date();
	return e = +e, t = +t, function(r) {
		return n.setTime(e * (1 - r) + t * r), n;
	};
}
//#endregion
//#region node_modules/d3-interpolate/src/number.js
function W(e, t) {
	return e = +e, t = +t, function(n) {
		return e * (1 - n) + t * n;
	};
}
//#endregion
//#region node_modules/d3-interpolate/src/object.js
function yr(e, t) {
	var n = {}, r = {}, i;
	for (i in (typeof e != "object" || !e) && (e = {}), (typeof t != "object" || !t) && (t = {}), t) i in e ? n[i] = Tr(e[i], t[i]) : r[i] = t[i];
	return function(e) {
		for (i in n) r[i] = n[i](e);
		return r;
	};
}
//#endregion
//#region node_modules/d3-interpolate/src/string.js
var br = /[-+]?(?:\d+\.?\d*|\.?\d+)(?:[eE][-+]?\d+)?/g, xr = new RegExp(br.source, "g");
function Sr(e) {
	return function() {
		return e;
	};
}
function Cr(e) {
	return function(t) {
		return e(t) + "";
	};
}
function wr(e, t) {
	var n = br.lastIndex = xr.lastIndex = 0, r, i, a, o = -1, s = [], c = [];
	for (e += "", t += ""; (r = br.exec(e)) && (i = xr.exec(t));) (a = i.index) > n && (a = t.slice(n, a), s[o] ? s[o] += a : s[++o] = a), (r = r[0]) === (i = i[0]) ? s[o] ? s[o] += i : s[++o] = i : (s[++o] = null, c.push({
		i: o,
		x: W(r, i)
	})), n = xr.lastIndex;
	return n < t.length && (a = t.slice(n), s[o] ? s[o] += a : s[++o] = a), s.length < 2 ? c[0] ? Cr(c[0].x) : Sr(t) : (t = c.length, function(e) {
		for (var n = 0, r; n < t; ++n) s[(r = c[n]).i] = r.x(e);
		return s.join("");
	});
}
//#endregion
//#region node_modules/d3-interpolate/src/value.js
function Tr(e, t) {
	var n = typeof t, r;
	return t == null || n === "boolean" ? lr(t) : (n === "number" ? W : n === "string" ? (r = Dn(t)) ? (t = r, mr) : wr : t instanceof Dn ? mr : t instanceof Date ? vr : gr(t) ? hr : Array.isArray(t) ? _r : typeof t.valueOf != "function" && typeof t.toString != "function" || isNaN(t) ? yr : W)(e, t);
}
//#endregion
//#region node_modules/d3-interpolate/src/round.js
function Er(e, t) {
	return e = +e, t = +t, function(n) {
		return Math.round(e * (1 - n) + t * n);
	};
}
//#endregion
//#region node_modules/d3-interpolate/src/transform/decompose.js
var Dr = 180 / Math.PI, Or = {
	translateX: 0,
	translateY: 0,
	rotate: 0,
	skewX: 0,
	scaleX: 1,
	scaleY: 1
};
function kr(e, t, n, r, i, a) {
	var o, s, c;
	return (o = Math.sqrt(e * e + t * t)) && (e /= o, t /= o), (c = e * n + t * r) && (n -= e * c, r -= t * c), (s = Math.sqrt(n * n + r * r)) && (n /= s, r /= s, c /= s), e * r < t * n && (e = -e, t = -t, c = -c, o = -o), {
		translateX: i,
		translateY: a,
		rotate: Math.atan2(t, e) * Dr,
		skewX: Math.atan(c) * Dr,
		scaleX: o,
		scaleY: s
	};
}
//#endregion
//#region node_modules/d3-interpolate/src/transform/parse.js
var Ar;
function jr(e) {
	let t = new (typeof DOMMatrix == "function" ? DOMMatrix : WebKitCSSMatrix)(e + "");
	return t.isIdentity ? Or : kr(t.a, t.b, t.c, t.d, t.e, t.f);
}
function Mr(e) {
	return e == null || (Ar ||= document.createElementNS("http://www.w3.org/2000/svg", "g"), Ar.setAttribute("transform", e), !(e = Ar.transform.baseVal.consolidate())) ? Or : (e = e.matrix, kr(e.a, e.b, e.c, e.d, e.e, e.f));
}
//#endregion
//#region node_modules/d3-interpolate/src/transform/index.js
function Nr(e, t, n, r) {
	function i(e) {
		return e.length ? e.pop() + " " : "";
	}
	function a(e, r, i, a, o, s) {
		if (e !== i || r !== a) {
			var c = o.push("translate(", null, t, null, n);
			s.push({
				i: c - 4,
				x: W(e, i)
			}, {
				i: c - 2,
				x: W(r, a)
			});
		} else (i || a) && o.push("translate(" + i + t + a + n);
	}
	function o(e, t, n, a) {
		e === t ? t && n.push(i(n) + "rotate(" + t + r) : (e - t > 180 ? t += 360 : t - e > 180 && (e += 360), a.push({
			i: n.push(i(n) + "rotate(", null, r) - 2,
			x: W(e, t)
		}));
	}
	function s(e, t, n, a) {
		e === t ? t && n.push(i(n) + "skewX(" + t + r) : a.push({
			i: n.push(i(n) + "skewX(", null, r) - 2,
			x: W(e, t)
		});
	}
	function c(e, t, n, r, a, o) {
		if (e !== n || t !== r) {
			var s = a.push(i(a) + "scale(", null, ",", null, ")");
			o.push({
				i: s - 4,
				x: W(e, n)
			}, {
				i: s - 2,
				x: W(t, r)
			});
		} else (n !== 1 || r !== 1) && a.push(i(a) + "scale(" + n + "," + r + ")");
	}
	return function(t, n) {
		var r = [], i = [];
		return t = e(t), n = e(n), a(t.translateX, t.translateY, n.translateX, n.translateY, r, i), o(t.rotate, n.rotate, r, i), s(t.skewX, n.skewX, r, i), c(t.scaleX, t.scaleY, n.scaleX, n.scaleY, r, i), t = n = null, function(e) {
			for (var t = -1, n = i.length, a; ++t < n;) r[(a = i[t]).i] = a.x(e);
			return r.join("");
		};
	};
}
var Pr = Nr(jr, "px, ", "px)", "deg)"), Fr = Nr(Mr, ", ", ")", ")"), Ir = 1e-12;
function Lr(e) {
	return ((e = Math.exp(e)) + 1 / e) / 2;
}
function Rr(e) {
	return ((e = Math.exp(e)) - 1 / e) / 2;
}
function zr(e) {
	return ((e = Math.exp(2 * e)) - 1) / (e + 1);
}
var Br = (function e(t, n, r) {
	function i(e, i) {
		var a = e[0], o = e[1], s = e[2], c = i[0], l = i[1], u = i[2], d = c - a, f = l - o, p = d * d + f * f, m, h;
		if (p < Ir) h = Math.log(u / s) / t, m = function(e) {
			return [
				a + e * d,
				o + e * f,
				s * Math.exp(t * e * h)
			];
		};
		else {
			var g = Math.sqrt(p), _ = (u * u - s * s + r * p) / (2 * s * n * g), v = (u * u - s * s - r * p) / (2 * u * n * g), y = Math.log(Math.sqrt(_ * _ + 1) - _);
			h = (Math.log(Math.sqrt(v * v + 1) - v) - y) / t, m = function(e) {
				var r = e * h, i = Lr(y), c = s / (n * g) * (i * zr(t * r + y) - Rr(y));
				return [
					a + c * d,
					o + c * f,
					s * i / Lr(t * r + y)
				];
			};
		}
		return m.duration = h * 1e3 * t / Math.SQRT2, m;
	}
	return i.rho = function(t) {
		var n = Math.max(.001, +t), r = n * n;
		return e(n, r, r * r);
	}, i;
})(Math.SQRT2, 2, 4);
//#endregion
//#region node_modules/d3-interpolate/src/lab.js
function Vr(e, t) {
	var n = pr((e = tr(e)).l, (t = tr(t)).l), r = pr(e.a, t.a), i = pr(e.b, t.b), a = pr(e.opacity, t.opacity);
	return function(t) {
		return e.l = n(t), e.a = r(t), e.b = i(t), e.opacity = a(t), e + "";
	};
}
//#endregion
//#region node_modules/d3-timer/src/timer.js
var Hr = 0, Ur = 0, Wr = 0, Gr = 1e3, Kr, qr, Jr = 0, Yr = 0, Xr = 0, Zr = typeof performance == "object" && performance.now ? performance : Date, Qr = typeof window == "object" && window.requestAnimationFrame ? window.requestAnimationFrame.bind(window) : function(e) {
	setTimeout(e, 17);
};
function $r() {
	return Yr ||= (Qr(ei), Zr.now() + Xr);
}
function ei() {
	Yr = 0;
}
function ti() {
	this._call = this._time = this._next = null;
}
ti.prototype = ni.prototype = {
	constructor: ti,
	restart: function(e, t, n) {
		if (typeof e != "function") throw TypeError("callback is not a function");
		n = (n == null ? $r() : +n) + (t == null ? 0 : +t), !this._next && qr !== this && (qr ? qr._next = this : Kr = this, qr = this), this._call = e, this._time = n, si();
	},
	stop: function() {
		this._call && (this._call = null, this._time = Infinity, si());
	}
};
function ni(e, t, n) {
	var r = new ti();
	return r.restart(e, t, n), r;
}
function ri() {
	$r(), ++Hr;
	for (var e = Kr, t; e;) (t = Yr - e._time) >= 0 && e._call.call(void 0, t), e = e._next;
	--Hr;
}
function ii() {
	Yr = (Jr = Zr.now()) + Xr, Hr = Ur = 0;
	try {
		ri();
	} finally {
		Hr = 0, oi(), Yr = 0;
	}
}
function ai() {
	var e = Zr.now(), t = e - Jr;
	t > Gr && (Xr -= t, Jr = e);
}
function oi() {
	for (var e, t = Kr, n, r = Infinity; t;) t._call ? (r > t._time && (r = t._time), e = t, t = t._next) : (n = t._next, t._next = null, t = e ? e._next = n : Kr = n);
	qr = e, si(r);
}
function si(e) {
	Hr || (Ur &&= clearTimeout(Ur), e - Yr > 24 ? (e < Infinity && (Ur = setTimeout(ii, e - Zr.now() - Xr)), Wr &&= clearInterval(Wr)) : (Wr ||= (Jr = Zr.now(), setInterval(ai, Gr)), Hr = 1, Qr(ii)));
}
//#endregion
//#region node_modules/d3-timer/src/timeout.js
function ci(e, t, n) {
	var r = new ti();
	return t = t == null ? 0 : +t, r.restart((n) => {
		r.stop(), e(n + t);
	}, t, n), r;
}
//#endregion
//#region node_modules/d3-transition/src/transition/schedule.js
var li = ne("start", "end", "cancel", "interrupt"), ui = [];
function di(e, t, n, r, i, a) {
	var o = e.__transition;
	if (!o) e.__transition = {};
	else if (n in o) return;
	pi(e, n, {
		name: t,
		index: r,
		group: i,
		on: li,
		tween: ui,
		time: a.time,
		delay: a.delay,
		duration: a.duration,
		ease: a.ease,
		timer: null,
		state: 0
	});
}
function fi(e, t) {
	var n = K(e, t);
	if (n.state > 0) throw Error("too late; already scheduled");
	return n;
}
function G(e, t) {
	var n = K(e, t);
	if (n.state > 3) throw Error("too late; already running");
	return n;
}
function K(e, t) {
	var n = e.__transition;
	if (!n || !(n = n[t])) throw Error("transition not found");
	return n;
}
function pi(e, t, n) {
	var r = e.__transition, i;
	r[t] = n, n.timer = ni(a, 0, n.time);
	function a(e) {
		n.state = 1, n.timer.restart(o, n.delay, n.time), n.delay <= e && o(e - n.delay);
	}
	function o(a) {
		var l, u, d, f;
		if (n.state !== 1) return c();
		for (l in r) if (f = r[l], f.name === n.name) {
			if (f.state === 3) return ci(o);
			f.state === 4 ? (f.state = 6, f.timer.stop(), f.on.call("interrupt", e, e.__data__, f.index, f.group), delete r[l]) : +l < t && (f.state = 6, f.timer.stop(), f.on.call("cancel", e, e.__data__, f.index, f.group), delete r[l]);
		}
		if (ci(function() {
			n.state === 3 && (n.state = 4, n.timer.restart(s, n.delay, n.time), s(a));
		}), n.state = 2, n.on.call("start", e, e.__data__, n.index, n.group), n.state === 2) {
			for (n.state = 3, i = Array(d = n.tween.length), l = 0, u = -1; l < d; ++l) (f = n.tween[l].value.call(e, e.__data__, n.index, n.group)) && (i[++u] = f);
			i.length = u + 1;
		}
	}
	function s(t) {
		for (var r = t < n.duration ? n.ease.call(null, t / n.duration) : (n.timer.restart(c), n.state = 5, 1), a = -1, o = i.length; ++a < o;) i[a].call(e, r);
		n.state === 5 && (n.on.call("end", e, e.__data__, n.index, n.group), c());
	}
	function c() {
		for (var i in n.state = 6, n.timer.stop(), delete r[t], r) return;
		delete e.__transition;
	}
}
//#endregion
//#region node_modules/d3-transition/src/interrupt.js
function mi(e, t) {
	var n = e.__transition, r, i, a = !0, o;
	if (n) {
		for (o in t = t == null ? null : t + "", n) {
			if ((r = n[o]).name !== t) {
				a = !1;
				continue;
			}
			i = r.state > 2 && r.state < 5, r.state = 6, r.timer.stop(), r.on.call(i ? "interrupt" : "cancel", e, e.__data__, r.index, r.group), delete n[o];
		}
		a && delete e.__transition;
	}
}
//#endregion
//#region node_modules/d3-transition/src/selection/interrupt.js
function hi(e) {
	return this.each(function() {
		mi(this, e);
	});
}
//#endregion
//#region node_modules/d3-transition/src/transition/tween.js
function gi(e, t) {
	var n, r;
	return function() {
		var i = G(this, e), a = i.tween;
		if (a !== n) {
			r = n = a;
			for (var o = 0, s = r.length; o < s; ++o) if (r[o].name === t) {
				r = r.slice(), r.splice(o, 1);
				break;
			}
		}
		i.tween = r;
	};
}
function _i(e, t, n) {
	var r, i;
	if (typeof n != "function") throw Error();
	return function() {
		var a = G(this, e), o = a.tween;
		if (o !== r) {
			i = (r = o).slice();
			for (var s = {
				name: t,
				value: n
			}, c = 0, l = i.length; c < l; ++c) if (i[c].name === t) {
				i[c] = s;
				break;
			}
			c === l && i.push(s);
		}
		a.tween = i;
	};
}
function vi(e, t) {
	var n = this._id;
	if (e += "", arguments.length < 2) {
		for (var r = K(this.node(), n).tween, i = 0, a = r.length, o; i < a; ++i) if ((o = r[i]).name === e) return o.value;
		return null;
	}
	return this.each((t == null ? gi : _i)(n, e, t));
}
function yi(e, t, n) {
	var r = e._id;
	return e.each(function() {
		var e = G(this, r);
		(e.value ||= {})[t] = n.apply(this, arguments);
	}), function(e) {
		return K(e, r).value[t];
	};
}
//#endregion
//#region node_modules/d3-transition/src/transition/interpolate.js
function bi(e, t) {
	var n;
	return (typeof t == "number" ? W : t instanceof Dn ? mr : (n = Dn(t)) ? (t = n, mr) : wr)(e, t);
}
//#endregion
//#region node_modules/d3-transition/src/transition/attr.js
function xi(e) {
	return function() {
		this.removeAttribute(e);
	};
}
function Si(e) {
	return function() {
		this.removeAttributeNS(e.space, e.local);
	};
}
function Ci(e, t, n) {
	var r, i = n + "", a;
	return function() {
		var o = this.getAttribute(e);
		return o === i ? null : o === r ? a : a = t(r = o, n);
	};
}
function wi(e, t, n) {
	var r, i = n + "", a;
	return function() {
		var o = this.getAttributeNS(e.space, e.local);
		return o === i ? null : o === r ? a : a = t(r = o, n);
	};
}
function Ti(e, t, n) {
	var r, i, a;
	return function() {
		var o, s = n(this), c;
		return s == null ? void this.removeAttribute(e) : (o = this.getAttribute(e), c = s + "", o === c ? null : o === r && c === i ? a : (i = c, a = t(r = o, s)));
	};
}
function Ei(e, t, n) {
	var r, i, a;
	return function() {
		var o, s = n(this), c;
		return s == null ? void this.removeAttributeNS(e.space, e.local) : (o = this.getAttributeNS(e.space, e.local), c = s + "", o === c ? null : o === r && c === i ? a : (i = c, a = t(r = o, s)));
	};
}
function Di(e, t) {
	var n = I(e), r = n === "transform" ? Fr : bi;
	return this.attrTween(e, typeof t == "function" ? (n.local ? Ei : Ti)(n, r, yi(this, "attr." + e, t)) : t == null ? (n.local ? Si : xi)(n) : (n.local ? wi : Ci)(n, r, t));
}
//#endregion
//#region node_modules/d3-transition/src/transition/attrTween.js
function Oi(e, t) {
	return function(n) {
		this.setAttribute(e, t.call(this, n));
	};
}
function ki(e, t) {
	return function(n) {
		this.setAttributeNS(e.space, e.local, t.call(this, n));
	};
}
function Ai(e, t) {
	var n, r;
	function i() {
		var i = t.apply(this, arguments);
		return i !== r && (n = (r = i) && ki(e, i)), n;
	}
	return i._value = t, i;
}
function ji(e, t) {
	var n, r;
	function i() {
		var i = t.apply(this, arguments);
		return i !== r && (n = (r = i) && Oi(e, i)), n;
	}
	return i._value = t, i;
}
function Mi(e, t) {
	var n = "attr." + e;
	if (arguments.length < 2) return (n = this.tween(n)) && n._value;
	if (t == null) return this.tween(n, null);
	if (typeof t != "function") throw Error();
	var r = I(e);
	return this.tween(n, (r.local ? Ai : ji)(r, t));
}
//#endregion
//#region node_modules/d3-transition/src/transition/delay.js
function Ni(e, t) {
	return function() {
		fi(this, e).delay = +t.apply(this, arguments);
	};
}
function Pi(e, t) {
	return t = +t, function() {
		fi(this, e).delay = t;
	};
}
function Fi(e) {
	var t = this._id;
	return arguments.length ? this.each((typeof e == "function" ? Ni : Pi)(t, e)) : K(this.node(), t).delay;
}
//#endregion
//#region node_modules/d3-transition/src/transition/duration.js
function Ii(e, t) {
	return function() {
		G(this, e).duration = +t.apply(this, arguments);
	};
}
function Li(e, t) {
	return t = +t, function() {
		G(this, e).duration = t;
	};
}
function Ri(e) {
	var t = this._id;
	return arguments.length ? this.each((typeof e == "function" ? Ii : Li)(t, e)) : K(this.node(), t).duration;
}
//#endregion
//#region node_modules/d3-transition/src/transition/ease.js
function zi(e, t) {
	if (typeof t != "function") throw Error();
	return function() {
		G(this, e).ease = t;
	};
}
function Bi(e) {
	var t = this._id;
	return arguments.length ? this.each(zi(t, e)) : K(this.node(), t).ease;
}
//#endregion
//#region node_modules/d3-transition/src/transition/easeVarying.js
function Vi(e, t) {
	return function() {
		var n = t.apply(this, arguments);
		if (typeof n != "function") throw Error();
		G(this, e).ease = n;
	};
}
function Hi(e) {
	if (typeof e != "function") throw Error();
	return this.each(Vi(this._id, e));
}
//#endregion
//#region node_modules/d3-transition/src/transition/filter.js
function Ui(e) {
	typeof e != "function" && (e = _e(e));
	for (var t = this._groups, n = t.length, r = Array(n), i = 0; i < n; ++i) for (var a = t[i], o = a.length, s = r[i] = [], c, l = 0; l < o; ++l) (c = a[l]) && e.call(c, c.__data__, l, a) && s.push(c);
	return new q(r, this._parents, this._name, this._id);
}
//#endregion
//#region node_modules/d3-transition/src/transition/merge.js
function Wi(e) {
	if (e._id !== this._id) throw Error();
	for (var t = this._groups, n = e._groups, r = t.length, i = n.length, a = Math.min(r, i), o = Array(r), s = 0; s < a; ++s) for (var c = t[s], l = n[s], u = c.length, d = o[s] = Array(u), f, p = 0; p < u; ++p) (f = c[p] || l[p]) && (d[p] = f);
	for (; s < r; ++s) o[s] = t[s];
	return new q(o, this._parents, this._name, this._id);
}
//#endregion
//#region node_modules/d3-transition/src/transition/on.js
function Gi(e) {
	return (e + "").trim().split(/^|\s+/).every(function(e) {
		var t = e.indexOf(".");
		return t >= 0 && (e = e.slice(0, t)), !e || e === "start";
	});
}
function Ki(e, t, n) {
	var r, i, a = Gi(t) ? fi : G;
	return function() {
		var o = a(this, e), s = o.on;
		s !== r && (i = (r = s).copy()).on(t, n), o.on = i;
	};
}
function qi(e, t) {
	var n = this._id;
	return arguments.length < 2 ? K(this.node(), n).on.on(e) : this.each(Ki(n, e, t));
}
//#endregion
//#region node_modules/d3-transition/src/transition/remove.js
function Ji(e) {
	return function() {
		var t = this.parentNode;
		for (var n in this.__transition) if (+n !== e) return;
		t && t.removeChild(this);
	};
}
function Yi() {
	return this.on("end.remove", Ji(this._id));
}
//#endregion
//#region node_modules/d3-transition/src/transition/select.js
function Xi(e) {
	var t = this._name, n = this._id;
	typeof e != "function" && (e = ue(e));
	for (var r = this._groups, i = r.length, a = Array(i), o = 0; o < i; ++o) for (var s = r[o], c = s.length, l = a[o] = Array(c), u, d, f = 0; f < c; ++f) (u = s[f]) && (d = e.call(u, u.__data__, f, s)) && ("__data__" in u && (d.__data__ = u.__data__), l[f] = d, di(l[f], t, n, f, l, K(u, n)));
	return new q(a, this._parents, t, n);
}
//#endregion
//#region node_modules/d3-transition/src/transition/selectAll.js
function Zi(e) {
	var t = this._name, n = this._id;
	typeof e != "function" && (e = me(e));
	for (var r = this._groups, i = r.length, a = [], o = [], s = 0; s < i; ++s) for (var c = r[s], l = c.length, u, d = 0; d < l; ++d) if (u = c[d]) {
		for (var f = e.call(u, u.__data__, d, c), p, m = K(u, n), h = 0, g = f.length; h < g; ++h) (p = f[h]) && di(p, t, n, h, f, m);
		a.push(f), o.push(u);
	}
	return new q(a, o, t, n);
}
//#endregion
//#region node_modules/d3-transition/src/transition/selection.js
var Qi = $t.prototype.constructor;
function $i() {
	return new Qi(this._groups, this._parents);
}
//#endregion
//#region node_modules/d3-transition/src/transition/style.js
function ea(e, t) {
	var n, r, i;
	return function() {
		var a = st(this, e), o = (this.style.removeProperty(e), st(this, e));
		return a === o ? null : a === n && o === r ? i : i = t(n = a, r = o);
	};
}
function ta(e) {
	return function() {
		this.style.removeProperty(e);
	};
}
function na(e, t, n) {
	var r, i = n + "", a;
	return function() {
		var o = st(this, e);
		return o === i ? null : o === r ? a : a = t(r = o, n);
	};
}
function ra(e, t, n) {
	var r, i, a;
	return function() {
		var o = st(this, e), s = n(this), c = s + "";
		return s ?? (c = s = (this.style.removeProperty(e), st(this, e))), o === c ? null : o === r && c === i ? a : (i = c, a = t(r = o, s));
	};
}
function ia(e, t) {
	var n, r, i, a = "style." + t, o = "end." + a, s;
	return function() {
		var c = G(this, e), l = c.on, u = c.value[a] == null ? s ||= ta(t) : void 0;
		(l !== n || i !== u) && (r = (n = l).copy()).on(o, i = u), c.on = r;
	};
}
function aa(e, t, n) {
	var r = (e += "") == "transform" ? Pr : bi;
	return t == null ? this.styleTween(e, ea(e, r)).on("end.style." + e, ta(e)) : typeof t == "function" ? this.styleTween(e, ra(e, r, yi(this, "style." + e, t))).each(ia(this._id, e)) : this.styleTween(e, na(e, r, t), n).on("end.style." + e, null);
}
//#endregion
//#region node_modules/d3-transition/src/transition/styleTween.js
function oa(e, t, n) {
	return function(r) {
		this.style.setProperty(e, t.call(this, r), n);
	};
}
function sa(e, t, n) {
	var r, i;
	function a() {
		var a = t.apply(this, arguments);
		return a !== i && (r = (i = a) && oa(e, a, n)), r;
	}
	return a._value = t, a;
}
function ca(e, t, n) {
	var r = "style." + (e += "");
	if (arguments.length < 2) return (r = this.tween(r)) && r._value;
	if (t == null) return this.tween(r, null);
	if (typeof t != "function") throw Error();
	return this.tween(r, sa(e, t, n ?? ""));
}
//#endregion
//#region node_modules/d3-transition/src/transition/text.js
function la(e) {
	return function() {
		this.textContent = e;
	};
}
function ua(e) {
	return function() {
		var t = e(this);
		this.textContent = t ?? "";
	};
}
function da(e) {
	return this.tween("text", typeof e == "function" ? ua(yi(this, "text", e)) : la(e == null ? "" : e + ""));
}
//#endregion
//#region node_modules/d3-transition/src/transition/textTween.js
function fa(e) {
	return function(t) {
		this.textContent = e.call(this, t);
	};
}
function pa(e) {
	var t, n;
	function r() {
		var r = e.apply(this, arguments);
		return r !== n && (t = (n = r) && fa(r)), t;
	}
	return r._value = e, r;
}
function ma(e) {
	var t = "text";
	if (arguments.length < 1) return (t = this.tween(t)) && t._value;
	if (e == null) return this.tween(t, null);
	if (typeof e != "function") throw Error();
	return this.tween(t, pa(e));
}
//#endregion
//#region node_modules/d3-transition/src/transition/transition.js
function ha() {
	for (var e = this._name, t = this._id, n = va(), r = this._groups, i = r.length, a = 0; a < i; ++a) for (var o = r[a], s = o.length, c, l = 0; l < s; ++l) if (c = o[l]) {
		var u = K(c, t);
		di(c, e, n, l, o, {
			time: u.time + u.delay + u.duration,
			delay: 0,
			duration: u.duration,
			ease: u.ease
		});
	}
	return new q(r, this._parents, e, n);
}
//#endregion
//#region node_modules/d3-transition/src/transition/end.js
function ga() {
	var e, t, n = this, r = n._id, i = n.size();
	return new Promise(function(a, o) {
		var s = { value: o }, c = { value: function() {
			--i === 0 && a();
		} };
		n.each(function() {
			var n = G(this, r), i = n.on;
			i !== e && (t = (e = i).copy(), t._.cancel.push(s), t._.interrupt.push(s), t._.end.push(c)), n.on = t;
		}), i === 0 && a();
	});
}
//#endregion
//#region node_modules/d3-transition/src/transition/index.js
var _a = 0;
function q(e, t, n, r) {
	this._groups = e, this._parents = t, this._name = n, this._id = r;
}
function va() {
	return ++_a;
}
var J = $t.prototype;
q.prototype = {
	constructor: q,
	select: Xi,
	selectAll: Zi,
	selectChild: J.selectChild,
	selectChildren: J.selectChildren,
	filter: Ui,
	merge: Wi,
	selection: $i,
	transition: ha,
	call: J.call,
	nodes: J.nodes,
	node: J.node,
	size: J.size,
	empty: J.empty,
	each: J.each,
	on: qi,
	attr: Di,
	attrTween: Mi,
	style: aa,
	styleTween: ca,
	text: da,
	textTween: ma,
	remove: Yi,
	tween: vi,
	delay: Fi,
	duration: Ri,
	ease: Bi,
	easeVarying: Hi,
	end: ga,
	[Symbol.iterator]: J[Symbol.iterator]
};
//#endregion
//#region node_modules/d3-ease/src/cubic.js
function ya(e) {
	return ((e *= 2) <= 1 ? e * e * e : (e -= 2) * e * e + 2) / 2;
}
//#endregion
//#region node_modules/d3-transition/src/selection/transition.js
var ba = {
	time: null,
	delay: 0,
	duration: 250,
	ease: ya
};
function xa(e, t) {
	for (var n; !(n = e.__transition) || !(n = n[t]);) if (!(e = e.parentNode)) throw Error(`transition ${t} not found`);
	return n;
}
function Sa(e) {
	var t, n;
	e instanceof q ? (t = e._id, e = e._name) : (t = va(), (n = ba).time = $r(), e = e == null ? null : e + "");
	for (var r = this._groups, i = r.length, a = 0; a < i; ++a) for (var o = r[a], s = o.length, c, l = 0; l < s; ++l) (c = o[l]) && di(c, e, t, l, o, n || xa(c, t));
	return new q(r, this._parents, e, t);
}
$t.prototype.interrupt = hi, $t.prototype.transition = Sa;
//#endregion
//#region node_modules/d3-brush/src/brush.js
var { abs: Ca, max: wa, min: Ta } = Math;
["w", "e"].map(Ea), ["n", "s"].map(Ea), [
	"n",
	"w",
	"e",
	"s",
	"nw",
	"ne",
	"sw",
	"se"
].map(Ea);
function Ea(e) {
	return { type: e };
}
//#endregion
//#region node_modules/d3-path/src/path.js
var Da = Math.PI, Oa = 2 * Da, ka = 1e-6, Aa = Oa - ka;
function ja(e) {
	this._ += e[0];
	for (let t = 1, n = e.length; t < n; ++t) this._ += arguments[t] + e[t];
}
function Ma(e) {
	let t = Math.floor(e);
	if (!(t >= 0)) throw Error(`invalid digits: ${e}`);
	if (t > 15) return ja;
	let n = 10 ** t;
	return function(e) {
		this._ += e[0];
		for (let t = 1, r = e.length; t < r; ++t) this._ += Math.round(arguments[t] * n) / n + e[t];
	};
}
var Na = class {
	constructor(e) {
		this._x0 = this._y0 = this._x1 = this._y1 = null, this._ = "", this._append = e == null ? ja : Ma(e);
	}
	moveTo(e, t) {
		this._append`M${this._x0 = this._x1 = +e},${this._y0 = this._y1 = +t}`;
	}
	closePath() {
		this._x1 !== null && (this._x1 = this._x0, this._y1 = this._y0, this._append`Z`);
	}
	lineTo(e, t) {
		this._append`L${this._x1 = +e},${this._y1 = +t}`;
	}
	quadraticCurveTo(e, t, n, r) {
		this._append`Q${+e},${+t},${this._x1 = +n},${this._y1 = +r}`;
	}
	bezierCurveTo(e, t, n, r, i, a) {
		this._append`C${+e},${+t},${+n},${+r},${this._x1 = +i},${this._y1 = +a}`;
	}
	arcTo(e, t, n, r, i) {
		if (e = +e, t = +t, n = +n, r = +r, i = +i, i < 0) throw Error(`negative radius: ${i}`);
		let a = this._x1, o = this._y1, s = n - e, c = r - t, l = a - e, u = o - t, d = l * l + u * u;
		if (this._x1 === null) this._append`M${this._x1 = e},${this._y1 = t}`;
		else if (d > ka) {
			if (!(Math.abs(u * s - c * l) > ka) || !i) this._append`L${this._x1 = e},${this._y1 = t}`;
			else {
				let f = n - a, p = r - o, m = s * s + c * c, h = f * f + p * p, g = Math.sqrt(m), _ = Math.sqrt(d), v = i * Math.tan((Da - Math.acos((m + d - h) / (2 * g * _))) / 2), y = v / _, b = v / g;
				Math.abs(y - 1) > ka && this._append`L${e + y * l},${t + y * u}`, this._append`A${i},${i},0,0,${+(u * f > l * p)},${this._x1 = e + b * s},${this._y1 = t + b * c}`;
			}
		}
	}
	arc(e, t, n, r, i, a) {
		if (e = +e, t = +t, n = +n, a = !!a, n < 0) throw Error(`negative radius: ${n}`);
		let o = n * Math.cos(r), s = n * Math.sin(r), c = e + o, l = t + s, u = 1 ^ a, d = a ? r - i : i - r;
		this._x1 === null ? this._append`M${c},${l}` : (Math.abs(this._x1 - c) > ka || Math.abs(this._y1 - l) > ka) && this._append`L${c},${l}`, n && (d < 0 && (d = d % Oa + Oa), d > Aa ? this._append`A${n},${n},0,1,${u},${e - o},${t - s}A${n},${n},0,1,${u},${this._x1 = c},${this._y1 = l}` : d > ka && this._append`A${n},${n},0,${+(d >= Da)},${u},${this._x1 = e + n * Math.cos(i)},${this._y1 = t + n * Math.sin(i)}`);
	}
	rect(e, t, n, r) {
		this._append`M${this._x0 = this._x1 = +e},${this._y0 = this._y1 = +t}h${n = +n}v${+r}h${-n}Z`;
	}
	toString() {
		return this._;
	}
};
Na.prototype;
//#endregion
//#region node_modules/d3-format/src/formatDecimal.js
function Pa(e) {
	return Math.abs(e = Math.round(e)) >= 1e21 ? e.toLocaleString("en").replace(/,/g, "") : e.toString(10);
}
function Fa(e, t) {
	if (!isFinite(e) || e === 0) return null;
	var n = (e = t ? e.toExponential(t - 1) : e.toExponential()).indexOf("e"), r = e.slice(0, n);
	return [r.length > 1 ? r[0] + r.slice(2) : r, +e.slice(n + 1)];
}
//#endregion
//#region node_modules/d3-format/src/exponent.js
function Ia(e) {
	return e = Fa(Math.abs(e)), e ? e[1] : NaN;
}
//#endregion
//#region node_modules/d3-format/src/formatGroup.js
function La(e, t) {
	return function(n, r) {
		for (var i = n.length, a = [], o = 0, s = e[0], c = 0; i > 0 && s > 0 && (c + s + 1 > r && (s = Math.max(1, r - c)), a.push(n.substring(i -= s, i + s)), !((c += s + 1) > r));) s = e[o = (o + 1) % e.length];
		return a.reverse().join(t);
	};
}
//#endregion
//#region node_modules/d3-format/src/formatNumerals.js
function Ra(e) {
	return function(t) {
		return t.replace(/[0-9]/g, function(t) {
			return e[+t];
		});
	};
}
//#endregion
//#region node_modules/d3-format/src/formatSpecifier.js
var za = /^(?:(.)?([<>=^]))?([+\-( ])?([$#])?(0)?(\d+)?(,)?(\.\d+)?(~)?([a-z%])?$/i;
function Ba(e) {
	if (!(t = za.exec(e))) throw Error("invalid format: " + e);
	var t;
	return new Va({
		fill: t[1],
		align: t[2],
		sign: t[3],
		symbol: t[4],
		zero: t[5],
		width: t[6],
		comma: t[7],
		precision: t[8] && t[8].slice(1),
		trim: t[9],
		type: t[10]
	});
}
Ba.prototype = Va.prototype;
function Va(e) {
	this.fill = e.fill === void 0 ? " " : e.fill + "", this.align = e.align === void 0 ? ">" : e.align + "", this.sign = e.sign === void 0 ? "-" : e.sign + "", this.symbol = e.symbol === void 0 ? "" : e.symbol + "", this.zero = !!e.zero, this.width = e.width === void 0 ? void 0 : +e.width, this.comma = !!e.comma, this.precision = e.precision === void 0 ? void 0 : +e.precision, this.trim = !!e.trim, this.type = e.type === void 0 ? "" : e.type + "";
}
Va.prototype.toString = function() {
	return this.fill + this.align + this.sign + this.symbol + (this.zero ? "0" : "") + (this.width === void 0 ? "" : Math.max(1, this.width | 0)) + (this.comma ? "," : "") + (this.precision === void 0 ? "" : "." + Math.max(0, this.precision | 0)) + (this.trim ? "~" : "") + this.type;
};
//#endregion
//#region node_modules/d3-format/src/formatTrim.js
function Ha(e) {
	out: for (var t = e.length, n = 1, r = -1, i; n < t; ++n) switch (e[n]) {
		case ".":
			r = i = n;
			break;
		case "0":
			r === 0 && (r = n), i = n;
			break;
		default:
			if (!+e[n]) break out;
			r > 0 && (r = 0);
	}
	return r > 0 ? e.slice(0, r) + e.slice(i + 1) : e;
}
//#endregion
//#region node_modules/d3-format/src/formatPrefixAuto.js
var Ua;
function Wa(e, t) {
	var n = Fa(e, t);
	if (!n) return Ua = void 0, e.toPrecision(t);
	var r = n[0], i = n[1], a = i - (Ua = Math.max(-8, Math.min(8, Math.floor(i / 3))) * 3) + 1, o = r.length;
	return a === o ? r : a > o ? r + Array(a - o + 1).join("0") : a > 0 ? r.slice(0, a) + "." + r.slice(a) : "0." + Array(1 - a).join("0") + Fa(e, Math.max(0, t + a - 1))[0];
}
//#endregion
//#region node_modules/d3-format/src/formatRounded.js
function Ga(e, t) {
	var n = Fa(e, t);
	if (!n) return e + "";
	var r = n[0], i = n[1];
	return i < 0 ? "0." + Array(-i).join("0") + r : r.length > i + 1 ? r.slice(0, i + 1) + "." + r.slice(i + 1) : r + Array(i - r.length + 2).join("0");
}
//#endregion
//#region node_modules/d3-format/src/formatTypes.js
var Ka = {
	"%": (e, t) => (e * 100).toFixed(t),
	b: (e) => Math.round(e).toString(2),
	c: (e) => e + "",
	d: Pa,
	e: (e, t) => e.toExponential(t),
	f: (e, t) => e.toFixed(t),
	g: (e, t) => e.toPrecision(t),
	o: (e) => Math.round(e).toString(8),
	p: (e, t) => Ga(e * 100, t),
	r: Ga,
	s: Wa,
	X: (e) => Math.round(e).toString(16).toUpperCase(),
	x: (e) => Math.round(e).toString(16)
};
//#endregion
//#region node_modules/d3-format/src/identity.js
function qa(e) {
	return e;
}
//#endregion
//#region node_modules/d3-format/src/locale.js
var Ja = Array.prototype.map, Ya = [
	"y",
	"z",
	"a",
	"f",
	"p",
	"n",
	"µ",
	"m",
	"",
	"k",
	"M",
	"G",
	"T",
	"P",
	"E",
	"Z",
	"Y"
];
function Xa(e) {
	var t = e.grouping === void 0 || e.thousands === void 0 ? qa : La(Ja.call(e.grouping, Number), e.thousands + ""), n = e.currency === void 0 ? "" : e.currency[0] + "", r = e.currency === void 0 ? "" : e.currency[1] + "", i = e.decimal === void 0 ? "." : e.decimal + "", a = e.numerals === void 0 ? qa : Ra(Ja.call(e.numerals, String)), o = e.percent === void 0 ? "%" : e.percent + "", s = e.minus === void 0 ? "−" : e.minus + "", c = e.nan === void 0 ? "NaN" : e.nan + "";
	function l(e, l) {
		e = Ba(e);
		var u = e.fill, d = e.align, f = e.sign, p = e.symbol, m = e.zero, h = e.width, g = e.comma, _ = e.precision, v = e.trim, y = e.type;
		y === "n" ? (g = !0, y = "g") : Ka[y] || (_ === void 0 && (_ = 12), v = !0, y = "g"), (m || u === "0" && d === "=") && (m = !0, u = "0", d = "=");
		var b = (l && l.prefix !== void 0 ? l.prefix : "") + (p === "$" ? n : p === "#" && /[boxX]/.test(y) ? "0" + y.toLowerCase() : ""), x = (p === "$" ? r : /[%p]/.test(y) ? o : "") + (l && l.suffix !== void 0 ? l.suffix : ""), S = Ka[y], C = /[defgprs%]/.test(y);
		_ = _ === void 0 ? 6 : /[gprs]/.test(y) ? Math.max(1, Math.min(21, _)) : Math.max(0, Math.min(20, _));
		function w(e) {
			var n = b, r = x, o, l, p;
			if (y === "c") r = S(e) + r, e = "";
			else {
				e = +e;
				var w = e < 0 || 1 / e < 0;
				if (e = isNaN(e) ? c : S(Math.abs(e), _), v && (e = Ha(e)), w && +e == 0 && f !== "+" && (w = !1), n = (w ? f === "(" ? f : s : f === "-" || f === "(" ? "" : f) + n, r = (y === "s" && !isNaN(e) && Ua !== void 0 ? Ya[8 + Ua / 3] : "") + r + (w && f === "(" ? ")" : ""), C) {
					for (o = -1, l = e.length; ++o < l;) if (p = e.charCodeAt(o), 48 > p || p > 57) {
						r = (p === 46 ? i + e.slice(o + 1) : e.slice(o)) + r, e = e.slice(0, o);
						break;
					}
				}
			}
			g && !m && (e = t(e, Infinity));
			var T = n.length + e.length + r.length, E = T < h ? Array(h - T + 1).join(u) : "";
			switch (g && m && (e = t(E + e, E.length ? h - r.length : Infinity), E = ""), d) {
				case "<":
					e = n + e + r + E;
					break;
				case "=":
					e = n + E + e + r;
					break;
				case "^":
					e = E.slice(0, T = E.length >> 1) + n + e + r + E.slice(T);
					break;
				default: e = E + n + e + r;
			}
			return a(e);
		}
		return w.toString = function() {
			return e + "";
		}, w;
	}
	function u(e, t) {
		var n = Math.max(-8, Math.min(8, Math.floor(Ia(t) / 3))) * 3, r = 10 ** -n, i = l((e = Ba(e), e.type = "f", e), { suffix: Ya[8 + n / 3] });
		return function(e) {
			return i(r * e);
		};
	}
	return {
		format: l,
		formatPrefix: u
	};
}
//#endregion
//#region node_modules/d3-format/src/defaultLocale.js
var Za, Qa, $a;
eo({
	thousands: ",",
	grouping: [3],
	currency: ["$", ""]
});
function eo(e) {
	return Za = Xa(e), Qa = Za.format, $a = Za.formatPrefix, Za;
}
//#endregion
//#region node_modules/d3-format/src/precisionFixed.js
function to(e) {
	return Math.max(0, -Ia(Math.abs(e)));
}
//#endregion
//#region node_modules/d3-format/src/precisionPrefix.js
function no(e, t) {
	return Math.max(0, Math.max(-8, Math.min(8, Math.floor(Ia(t) / 3))) * 3 - Ia(Math.abs(e)));
}
//#endregion
//#region node_modules/d3-format/src/precisionRound.js
function ro(e, t) {
	return e = Math.abs(e), t = Math.abs(t) - e, Math.max(0, Ia(t) - Ia(e)) + 1;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/count.js
function io(e) {
	var t = 0, n = e.children, r = n && n.length;
	if (!r) t = 1;
	else for (; --r >= 0;) t += n[r].value;
	e.value = t;
}
function ao() {
	return this.eachAfter(io);
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/each.js
function oo(e, t) {
	let n = -1;
	for (let r of this) e.call(t, r, ++n, this);
	return this;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/eachBefore.js
function so(e, t) {
	for (var n = this, r = [n], i, a, o = -1; n = r.pop();) if (e.call(t, n, ++o, this), i = n.children) for (a = i.length - 1; a >= 0; --a) r.push(i[a]);
	return this;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/eachAfter.js
function co(e, t) {
	for (var n = this, r = [n], i = [], a, o, s, c = -1; n = r.pop();) if (i.push(n), a = n.children) for (o = 0, s = a.length; o < s; ++o) r.push(a[o]);
	for (; n = i.pop();) e.call(t, n, ++c, this);
	return this;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/find.js
function lo(e, t) {
	let n = -1;
	for (let r of this) if (e.call(t, r, ++n, this)) return r;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/sum.js
function uo(e) {
	return this.eachAfter(function(t) {
		for (var n = +e(t.data) || 0, r = t.children, i = r && r.length; --i >= 0;) n += r[i].value;
		t.value = n;
	});
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/sort.js
function fo(e) {
	return this.eachBefore(function(t) {
		t.children && t.children.sort(e);
	});
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/path.js
function po(e) {
	for (var t = this, n = mo(t, e), r = [t]; t !== n;) t = t.parent, r.push(t);
	for (var i = r.length; e !== n;) r.splice(i, 0, e), e = e.parent;
	return r;
}
function mo(e, t) {
	if (e === t) return e;
	var n = e.ancestors(), r = t.ancestors(), i = null;
	for (e = n.pop(), t = r.pop(); e === t;) i = e, e = n.pop(), t = r.pop();
	return i;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/ancestors.js
function ho() {
	for (var e = this, t = [e]; e = e.parent;) t.push(e);
	return t;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/descendants.js
function go() {
	return Array.from(this);
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/leaves.js
function _o() {
	var e = [];
	return this.eachBefore(function(t) {
		t.children || e.push(t);
	}), e;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/links.js
function vo() {
	var e = this, t = [];
	return e.each(function(n) {
		n !== e && t.push({
			source: n.parent,
			target: n
		});
	}), t;
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/iterator.js
function* yo() {
	var e = this, t, n = [e], r, i, a;
	do
		for (t = n.reverse(), n = []; e = t.pop();) if (yield e, r = e.children) for (i = 0, a = r.length; i < a; ++i) n.push(r[i]);
	while (n.length);
}
//#endregion
//#region node_modules/d3-hierarchy/src/hierarchy/index.js
function bo(e, t) {
	e instanceof Map ? (e = [void 0, e], t === void 0 && (t = Co)) : t === void 0 && (t = So);
	for (var n = new Eo(e), r, i = [n], a, o, s, c; r = i.pop();) if ((o = t(r.data)) && (c = (o = Array.from(o)).length)) for (r.children = o, s = c - 1; s >= 0; --s) i.push(a = o[s] = new Eo(o[s])), a.parent = r, a.depth = r.depth + 1;
	return n.eachBefore(To);
}
function xo() {
	return bo(this).eachBefore(wo);
}
function So(e) {
	return e.children;
}
function Co(e) {
	return Array.isArray(e) ? e[1] : null;
}
function wo(e) {
	e.data.value !== void 0 && (e.value = e.data.value), e.data = e.data.data;
}
function To(e) {
	var t = 0;
	do
		e.height = t;
	while ((e = e.parent) && e.height < ++t);
}
function Eo(e) {
	this.data = e, this.depth = this.height = 0, this.parent = null;
}
Eo.prototype = bo.prototype = {
	constructor: Eo,
	count: ao,
	each: oo,
	eachAfter: co,
	eachBefore: so,
	find: lo,
	sum: uo,
	sort: fo,
	path: po,
	ancestors: ho,
	descendants: go,
	leaves: _o,
	links: vo,
	copy: xo,
	[Symbol.iterator]: yo
};
//#endregion
//#region node_modules/d3-hierarchy/src/accessors.js
function Do(e) {
	if (typeof e != "function") throw Error();
	return e;
}
//#endregion
//#region node_modules/d3-hierarchy/src/constant.js
function Oo() {
	return 0;
}
function ko(e) {
	return function() {
		return e;
	};
}
//#endregion
//#region node_modules/d3-hierarchy/src/treemap/round.js
function Ao(e) {
	e.x0 = Math.round(e.x0), e.y0 = Math.round(e.y0), e.x1 = Math.round(e.x1), e.y1 = Math.round(e.y1);
}
//#endregion
//#region node_modules/d3-hierarchy/src/treemap/dice.js
function jo(e, t, n, r, i) {
	for (var a = e.children, o, s = -1, c = a.length, l = e.value && (r - t) / e.value; ++s < c;) o = a[s], o.y0 = n, o.y1 = i, o.x0 = t, o.x1 = t += o.value * l;
}
//#endregion
//#region node_modules/d3-hierarchy/src/partition.js
function Mo() {
	var e = 1, t = 1, n = 0, r = !1;
	function i(i) {
		var o = i.height + 1;
		return i.x0 = i.y0 = n, i.x1 = e, i.y1 = t / o, i.eachBefore(a(t, o)), r && i.eachBefore(Ao), i;
	}
	function a(e, t) {
		return function(r) {
			r.children && jo(r, r.x0, e * (r.depth + 1) / t, r.x1, e * (r.depth + 2) / t);
			var i = r.x0, a = r.y0, o = r.x1 - n, s = r.y1 - n;
			o < i && (i = o = (i + o) / 2), s < a && (a = s = (a + s) / 2), r.x0 = i, r.y0 = a, r.x1 = o, r.y1 = s;
		};
	}
	return i.round = function(e) {
		return arguments.length ? (r = !!e, i) : r;
	}, i.size = function(n) {
		return arguments.length ? (e = +n[0], t = +n[1], i) : [e, t];
	}, i.padding = function(e) {
		return arguments.length ? (n = +e, i) : n;
	}, i;
}
//#endregion
//#region node_modules/d3-hierarchy/src/tree.js
function No(e, t) {
	return e.parent === t.parent ? 1 : 2;
}
function Po(e) {
	var t = e.children;
	return t ? t[0] : e.t;
}
function Fo(e) {
	var t = e.children;
	return t ? t[t.length - 1] : e.t;
}
function Io(e, t, n) {
	var r = n / (t.i - e.i);
	t.c -= r, t.s += n, e.c += r, t.z += n, t.m += n;
}
function Lo(e) {
	for (var t = 0, n = 0, r = e.children, i = r.length, a; --i >= 0;) a = r[i], a.z += t, a.m += t, t += a.s + (n += a.c);
}
function Ro(e, t, n) {
	return e.a.parent === t.parent ? e.a : n;
}
function zo(e, t) {
	this._ = e, this.parent = null, this.children = null, this.A = null, this.a = this, this.z = 0, this.m = 0, this.c = 0, this.s = 0, this.t = null, this.i = t;
}
zo.prototype = Object.create(Eo.prototype);
function Bo(e) {
	for (var t = new zo(e, 0), n, r = [t], i, a, o, s; n = r.pop();) if (a = n._.children) for (n.children = Array(s = a.length), o = s - 1; o >= 0; --o) r.push(i = n.children[o] = new zo(a[o], o)), i.parent = n;
	return (t.parent = new zo(null, 0)).children = [t], t;
}
function Vo() {
	var e = No, t = 1, n = 1, r = null;
	function i(i) {
		var s = Bo(i);
		if (s.eachAfter(a), s.parent.m = -s.z, s.eachBefore(o), r) i.eachBefore(c);
		else {
			var l = i, u = i, d = i;
			i.eachBefore(function(e) {
				e.x < l.x && (l = e), e.x > u.x && (u = e), e.depth > d.depth && (d = e);
			});
			var f = l === u ? 1 : e(l, u) / 2, p = f - l.x, m = t / (u.x + f + p), h = n / (d.depth || 1);
			i.eachBefore(function(e) {
				e.x = (e.x + p) * m, e.y = e.depth * h;
			});
		}
		return i;
	}
	function a(t) {
		var n = t.children, r = t.parent.children, i = t.i ? r[t.i - 1] : null;
		if (n) {
			Lo(t);
			var a = (n[0].z + n[n.length - 1].z) / 2;
			i ? (t.z = i.z + e(t._, i._), t.m = t.z - a) : t.z = a;
		} else i && (t.z = i.z + e(t._, i._));
		t.parent.A = s(t, i, t.parent.A || r[0]);
	}
	function o(e) {
		e._.x = e.z + e.parent.m, e.m += e.parent.m;
	}
	function s(t, n, r) {
		if (n) {
			for (var i = t, a = t, o = n, s = i.parent.children[0], c = i.m, l = a.m, u = o.m, d = s.m, f; o = Fo(o), i = Po(i), o && i;) s = Po(s), a = Fo(a), a.a = t, f = o.z + u - i.z - c + e(o._, i._), f > 0 && (Io(Ro(o, t, r), t, f), c += f, l += f), u += o.m, c += i.m, d += s.m, l += a.m;
			o && !Fo(a) && (a.t = o, a.m += u - l), i && !Po(s) && (s.t = i, s.m += c - d, r = t);
		}
		return r;
	}
	function c(e) {
		e.x *= t, e.y = e.depth * n;
	}
	return i.separation = function(t) {
		return arguments.length ? (e = t, i) : e;
	}, i.size = function(e) {
		return arguments.length ? (r = !1, t = +e[0], n = +e[1], i) : r ? null : [t, n];
	}, i.nodeSize = function(e) {
		return arguments.length ? (r = !0, t = +e[0], n = +e[1], i) : r ? [t, n] : null;
	}, i;
}
//#endregion
//#region node_modules/d3-hierarchy/src/treemap/slice.js
function Ho(e, t, n, r, i) {
	for (var a = e.children, o, s = -1, c = a.length, l = e.value && (i - n) / e.value; ++s < c;) o = a[s], o.x0 = t, o.x1 = r, o.y0 = n, o.y1 = n += o.value * l;
}
//#endregion
//#region node_modules/d3-hierarchy/src/treemap/squarify.js
var Uo = (1 + Math.sqrt(5)) / 2;
function Wo(e, t, n, r, i, a) {
	for (var o = [], s = t.children, c, l, u = 0, d = 0, f = s.length, p, m, h = t.value, g, _, v, y, b, x, S; u < f;) {
		p = i - n, m = a - r;
		do
			g = s[d++].value;
		while (!g && d < f);
		for (_ = v = g, x = Math.max(m / p, p / m) / (h * e), S = g * g * x, b = Math.max(v / S, S / _); d < f; ++d) {
			if (g += l = s[d].value, l < _ && (_ = l), l > v && (v = l), S = g * g * x, y = Math.max(v / S, S / _), y > b) {
				g -= l;
				break;
			}
			b = y;
		}
		o.push(c = {
			value: g,
			dice: p < m,
			children: s.slice(u, d)
		}), c.dice ? jo(c, n, r, i, h ? r += m * g / h : a) : Ho(c, n, r, h ? n += p * g / h : i, a), h -= g, u = d;
	}
	return o;
}
var Go = (function e(t) {
	function n(e, n, r, i, a) {
		Wo(t, e, n, r, i, a);
	}
	return n.ratio = function(t) {
		return e((t = +t) > 1 ? t : 1);
	}, n;
})(Uo);
//#endregion
//#region node_modules/d3-hierarchy/src/treemap/index.js
function Ko() {
	var e = Go, t = !1, n = 1, r = 1, i = [0], a = Oo, o = Oo, s = Oo, c = Oo, l = Oo;
	function u(e) {
		return e.x0 = e.y0 = 0, e.x1 = n, e.y1 = r, e.eachBefore(d), i = [0], t && e.eachBefore(Ao), e;
	}
	function d(t) {
		var n = i[t.depth], r = t.x0 + n, u = t.y0 + n, d = t.x1 - n, f = t.y1 - n;
		d < r && (r = d = (r + d) / 2), f < u && (u = f = (u + f) / 2), t.x0 = r, t.y0 = u, t.x1 = d, t.y1 = f, t.children && (n = i[t.depth + 1] = a(t) / 2, r += l(t) - n, u += o(t) - n, d -= s(t) - n, f -= c(t) - n, d < r && (r = d = (r + d) / 2), f < u && (u = f = (u + f) / 2), e(t, r, u, d, f));
	}
	return u.round = function(e) {
		return arguments.length ? (t = !!e, u) : t;
	}, u.size = function(e) {
		return arguments.length ? (n = +e[0], r = +e[1], u) : [n, r];
	}, u.tile = function(t) {
		return arguments.length ? (e = Do(t), u) : e;
	}, u.padding = function(e) {
		return arguments.length ? u.paddingInner(e).paddingOuter(e) : u.paddingInner();
	}, u.paddingInner = function(e) {
		return arguments.length ? (a = typeof e == "function" ? e : ko(+e), u) : a;
	}, u.paddingOuter = function(e) {
		return arguments.length ? u.paddingTop(e).paddingRight(e).paddingBottom(e).paddingLeft(e) : u.paddingTop();
	}, u.paddingTop = function(e) {
		return arguments.length ? (o = typeof e == "function" ? e : ko(+e), u) : o;
	}, u.paddingRight = function(e) {
		return arguments.length ? (s = typeof e == "function" ? e : ko(+e), u) : s;
	}, u.paddingBottom = function(e) {
		return arguments.length ? (c = typeof e == "function" ? e : ko(+e), u) : c;
	}, u.paddingLeft = function(e) {
		return arguments.length ? (l = typeof e == "function" ? e : ko(+e), u) : l;
	}, u;
}
//#endregion
//#region node_modules/d3-scale/src/init.js
function qo(e, t) {
	switch (arguments.length) {
		case 0: break;
		case 1:
			this.range(e);
			break;
		default: this.range(t).domain(e);
	}
	return this;
}
//#endregion
//#region node_modules/d3-scale/src/ordinal.js
var Jo = Symbol("implicit");
function Yo() {
	var e = new c(), t = [], n = [], r = Jo;
	function i(i) {
		let a = e.get(i);
		if (a === void 0) {
			if (r !== Jo) return r;
			e.set(i, a = t.push(i) - 1);
		}
		return n[a % n.length];
	}
	return i.domain = function(n) {
		if (!arguments.length) return t.slice();
		t = [], e = new c();
		for (let r of n) e.has(r) || e.set(r, t.push(r) - 1);
		return i;
	}, i.range = function(e) {
		return arguments.length ? (n = Array.from(e), i) : n.slice();
	}, i.unknown = function(e) {
		return arguments.length ? (r = e, i) : r;
	}, i.copy = function() {
		return Yo(t, n).unknown(r);
	}, qo.apply(i, arguments), i;
}
//#endregion
//#region node_modules/d3-scale/src/band.js
function Xo() {
	var e = Yo().unknown(void 0), t = e.domain, n = e.range, r = 0, i = 1, a, o, s = !1, c = 0, l = 0, u = .5;
	delete e.unknown;
	function d() {
		var e = t().length, d = i < r, f = d ? i : r, p = d ? r : i;
		a = (p - f) / Math.max(1, e - c + l * 2), s && (a = Math.floor(a)), f += (p - f - a * (e - c)) * u, o = a * (1 - c), s && (f = Math.round(f), o = Math.round(o));
		var m = x(e).map(function(e) {
			return f + a * e;
		});
		return n(d ? m.reverse() : m);
	}
	return e.domain = function(e) {
		return arguments.length ? (t(e), d()) : t();
	}, e.range = function(e) {
		return arguments.length ? ([r, i] = e, r = +r, i = +i, d()) : [r, i];
	}, e.rangeRound = function(e) {
		return [r, i] = e, r = +r, i = +i, s = !0, d();
	}, e.bandwidth = function() {
		return o;
	}, e.step = function() {
		return a;
	}, e.round = function(e) {
		return arguments.length ? (s = !!e, d()) : s;
	}, e.padding = function(e) {
		return arguments.length ? (c = Math.min(1, l = +e), d()) : c;
	}, e.paddingInner = function(e) {
		return arguments.length ? (c = Math.min(1, e), d()) : c;
	}, e.paddingOuter = function(e) {
		return arguments.length ? (l = +e, d()) : l;
	}, e.align = function(e) {
		return arguments.length ? (u = Math.max(0, Math.min(1, e)), d()) : u;
	}, e.copy = function() {
		return Xo(t(), [r, i]).round(s).paddingInner(c).paddingOuter(l).align(u);
	}, qo.apply(d(), arguments);
}
//#endregion
//#region node_modules/d3-scale/src/constant.js
function Zo(e) {
	return function() {
		return e;
	};
}
//#endregion
//#region node_modules/d3-scale/src/number.js
function Qo(e) {
	return +e;
}
//#endregion
//#region node_modules/d3-scale/src/continuous.js
var $o = [0, 1];
function es(e) {
	return e;
}
function ts(e, t) {
	return (t -= e = +e) ? function(n) {
		return (n - e) / t;
	} : Zo(isNaN(t) ? NaN : .5);
}
function ns(e, t) {
	var n;
	return e > t && (n = e, e = t, t = n), function(n) {
		return Math.max(e, Math.min(t, n));
	};
}
function rs(e, t, n) {
	var r = e[0], i = e[1], a = t[0], o = t[1];
	return i < r ? (r = ts(i, r), a = n(o, a)) : (r = ts(r, i), a = n(a, o)), function(e) {
		return a(r(e));
	};
}
function is(e, t, n) {
	var r = Math.min(e.length, t.length) - 1, i = Array(r), a = Array(r), o = -1;
	for (e[r] < e[0] && (e = e.slice().reverse(), t = t.slice().reverse()); ++o < r;) i[o] = ts(e[o], e[o + 1]), a[o] = n(t[o], t[o + 1]);
	return function(t) {
		var n = s(e, t, 1, r) - 1;
		return a[n](i[n](t));
	};
}
function as(e, t) {
	return t.domain(e.domain()).range(e.range()).interpolate(e.interpolate()).clamp(e.clamp()).unknown(e.unknown());
}
function os() {
	var e = $o, t = $o, n = Tr, r, i, a, o = es, s, c, l;
	function u() {
		var n = Math.min(e.length, t.length);
		return o !== es && (o = ns(e[0], e[n - 1])), s = n > 2 ? is : rs, c = l = null, d;
	}
	function d(i) {
		return i == null || isNaN(i = +i) ? a : (c ||= s(e.map(r), t, n))(r(o(i)));
	}
	return d.invert = function(n) {
		return o(i((l ||= s(t, e.map(r), W))(n)));
	}, d.domain = function(t) {
		return arguments.length ? (e = Array.from(t, Qo), u()) : e.slice();
	}, d.range = function(e) {
		return arguments.length ? (t = Array.from(e), u()) : t.slice();
	}, d.rangeRound = function(e) {
		return t = Array.from(e), n = Er, u();
	}, d.clamp = function(e) {
		return arguments.length ? (o = e ? !0 : es, u()) : o !== es;
	}, d.interpolate = function(e) {
		return arguments.length ? (n = e, u()) : n;
	}, d.unknown = function(e) {
		return arguments.length ? (a = e, d) : a;
	}, function(e, t) {
		return r = e, i = t, u();
	};
}
function ss() {
	return os()(es, es);
}
//#endregion
//#region node_modules/d3-scale/src/tickFormat.js
function cs(e, t, n, r) {
	var i = y(e, t, n), a;
	switch (r = Ba(r ?? ",f"), r.type) {
		case "s":
			var o = Math.max(Math.abs(e), Math.abs(t));
			return r.precision == null && !isNaN(a = no(i, o)) && (r.precision = a), $a(r, o);
		case "":
		case "e":
		case "g":
		case "p":
		case "r":
			r.precision == null && !isNaN(a = ro(i, Math.max(Math.abs(e), Math.abs(t)))) && (r.precision = a - (r.type === "e"));
			break;
		case "f":
		case "%": r.precision == null && !isNaN(a = to(i)) && (r.precision = a - (r.type === "%") * 2);
	}
	return Qa(r);
}
//#endregion
//#region node_modules/d3-scale/src/linear.js
function ls(e) {
	var t = e.domain;
	return e.ticks = function(e) {
		var n = t();
		return _(n[0], n[n.length - 1], e ?? 10);
	}, e.tickFormat = function(e, n) {
		var r = t();
		return cs(r[0], r[r.length - 1], e ?? 10, n);
	}, e.nice = function(n) {
		n ??= 10;
		var r = t(), i = 0, a = r.length - 1, o = r[i], s = r[a], c, l, u = 10;
		for (s < o && (l = o, o = s, s = l, l = i, i = a, a = l); u-- > 0;) {
			if (l = v(o, s, n), l === c) return r[i] = o, r[a] = s, t(r);
			if (l > 0) o = Math.floor(o / l) * l, s = Math.ceil(s / l) * l;
			else if (l < 0) o = Math.ceil(o * l) / l, s = Math.floor(s * l) / l;
			else break;
			c = l;
		}
		return e;
	}, e;
}
function Y() {
	var e = ss();
	return e.copy = function() {
		return as(e, Y());
	}, qo.apply(e, arguments), ls(e);
}
//#endregion
//#region node_modules/d3-scale/src/quantize.js
function us() {
	var e = 0, t = 1, n = 1, r = [.5], i = [0, 1], a;
	function o(e) {
		return e != null && e <= e ? i[s(r, e, 0, n)] : a;
	}
	function c() {
		var i = -1;
		for (r = Array(n); ++i < n;) r[i] = ((i + 1) * t - (i - n) * e) / (n + 1);
		return o;
	}
	return o.domain = function(n) {
		return arguments.length ? ([e, t] = n, e = +e, t = +t, c()) : [e, t];
	}, o.range = function(e) {
		return arguments.length ? (n = (i = Array.from(e)).length - 1, c()) : i.slice();
	}, o.invertExtent = function(a) {
		var o = i.indexOf(a);
		return o < 0 ? [NaN, NaN] : o < 1 ? [e, r[0]] : o >= n ? [r[n - 1], t] : [r[o - 1], r[o]];
	}, o.unknown = function(e) {
		return arguments.length && (a = e), o;
	}, o.thresholds = function() {
		return r.slice();
	}, o.copy = function() {
		return us().domain([e, t]).range(i).unknown(a);
	}, qo.apply(ls(o), arguments);
}
//#endregion
//#region node_modules/d3-scale-chromatic/src/colors.js
function ds(e) {
	for (var t = e.length / 6 | 0, n = Array(t), r = 0; r < t;) n[r] = "#" + e.slice(r * 6, ++r * 6);
	return n;
}
//#endregion
//#region node_modules/d3-scale-chromatic/src/categorical/category10.js
var fs = ds("1f77b4ff7f0e2ca02cd627289467bd8c564be377c27f7f7fbcbd2217becf");
//#endregion
//#region node_modules/d3-shape/src/constant.js
function X(e) {
	return function() {
		return e;
	};
}
//#endregion
//#region node_modules/d3-shape/src/math.js
var ps = Math.abs, Z = Math.atan2, ms = Math.cos, hs = Math.max, gs = Math.min, Q = Math.sin, _s = Math.sqrt, vs = Math.PI, ys = vs / 2, bs = 2 * vs;
function xs(e) {
	return e > 1 ? 0 : e < -1 ? vs : Math.acos(e);
}
function Ss(e) {
	return e >= 1 ? ys : e <= -1 ? -ys : Math.asin(e);
}
//#endregion
//#region node_modules/d3-shape/src/path.js
function Cs(e) {
	let t = 3;
	return e.digits = function(n) {
		if (!arguments.length) return t;
		if (n == null) t = null;
		else {
			let e = Math.floor(n);
			if (!(e >= 0)) throw RangeError(`invalid digits: ${n}`);
			t = e;
		}
		return e;
	}, () => new Na(t);
}
//#endregion
//#region node_modules/d3-shape/src/arc.js
function ws(e) {
	return e.innerRadius;
}
function Ts(e) {
	return e.outerRadius;
}
function Es(e) {
	return e.startAngle;
}
function Ds(e) {
	return e.endAngle;
}
function Os(e) {
	return e && e.padAngle;
}
function ks(e, t, n, r, i, a, o, s) {
	var c = n - e, l = r - t, u = o - i, d = s - a, f = d * c - u * l;
	if (!(f * f < 1e-12)) return f = (u * (t - a) - d * (e - i)) / f, [e + f * c, t + f * l];
}
function As(e, t, n, r, i, a, o) {
	var s = e - n, c = t - r, l = (o ? a : -a) / _s(s * s + c * c), u = l * c, d = -l * s, f = e + u, p = t + d, m = n + u, h = r + d, g = (f + m) / 2, _ = (p + h) / 2, v = m - f, y = h - p, b = v * v + y * y, x = i - a, S = f * h - m * p, C = (y < 0 ? -1 : 1) * _s(hs(0, x * x * b - S * S)), w = (S * y - v * C) / b, T = (-S * v - y * C) / b, E = (S * y + v * C) / b, D = (-S * v + y * C) / b, O = w - g, k = T - _, A = E - g, ee = D - _;
	return O * O + k * k > A * A + ee * ee && (w = E, T = D), {
		cx: w,
		cy: T,
		x01: -u,
		y01: -d,
		x11: w * (i / x - 1),
		y11: T * (i / x - 1)
	};
}
function js() {
	var e = ws, t = Ts, n = X(0), r = null, i = Es, a = Ds, o = Os, s = null, c = Cs(l);
	function l() {
		var l, u, d = +e.apply(this, arguments), f = +t.apply(this, arguments), p = i.apply(this, arguments) - ys, m = a.apply(this, arguments) - ys, h = ps(m - p), g = m > p;
		if (s ||= l = c(), f < d && (u = f, f = d, d = u), !(f > 1e-12)) s.moveTo(0, 0);
		else if (h > bs - 1e-12) s.moveTo(f * ms(p), f * Q(p)), s.arc(0, 0, f, p, m, !g), d > 1e-12 && (s.moveTo(d * ms(m), d * Q(m)), s.arc(0, 0, d, m, p, g));
		else {
			var _ = p, v = m, y = p, b = m, x = h, S = h, C = o.apply(this, arguments) / 2, w = C > 1e-12 && (r ? +r.apply(this, arguments) : _s(d * d + f * f)), T = gs(ps(f - d) / 2, +n.apply(this, arguments)), E = T, D = T, O, k;
			if (w > 1e-12) {
				var A = Ss(w / d * Q(C)), ee = Ss(w / f * Q(C));
				(x -= A * 2) > 1e-12 ? (A *= g ? 1 : -1, y += A, b -= A) : (x = 0, y = b = (p + m) / 2), (S -= ee * 2) > 1e-12 ? (ee *= g ? 1 : -1, _ += ee, v -= ee) : (S = 0, _ = v = (p + m) / 2);
			}
			var te = f * ms(_), j = f * Q(_), M = d * ms(b), N = d * Q(b);
			if (T > 1e-12) {
				var ne = f * ms(v), P = f * Q(v), re = d * ms(y), ie = d * Q(y), F;
				if (h < vs) {
					if (F = ks(te, j, re, ie, ne, P, M, N)) {
						var ae = te - F[0], I = j - F[1], oe = ne - F[0], se = P - F[1], ce = 1 / Q(xs((ae * oe + I * se) / (_s(ae * ae + I * I) * _s(oe * oe + se * se))) / 2), le = _s(F[0] * F[0] + F[1] * F[1]);
						E = gs(T, (d - le) / (ce - 1)), D = gs(T, (f - le) / (ce + 1));
					} else E = D = 0;
				}
			}
			S > 1e-12 ? D > 1e-12 ? (O = As(re, ie, te, j, f, D, g), k = As(ne, P, M, N, f, D, g), s.moveTo(O.cx + O.x01, O.cy + O.y01), D < T ? s.arc(O.cx, O.cy, D, Z(O.y01, O.x01), Z(k.y01, k.x01), !g) : (s.arc(O.cx, O.cy, D, Z(O.y01, O.x01), Z(O.y11, O.x11), !g), s.arc(0, 0, f, Z(O.cy + O.y11, O.cx + O.x11), Z(k.cy + k.y11, k.cx + k.x11), !g), s.arc(k.cx, k.cy, D, Z(k.y11, k.x11), Z(k.y01, k.x01), !g))) : (s.moveTo(te, j), s.arc(0, 0, f, _, v, !g)) : s.moveTo(te, j), !(d > 1e-12) || !(x > 1e-12) ? s.lineTo(M, N) : E > 1e-12 ? (O = As(M, N, ne, P, d, -E, g), k = As(te, j, re, ie, d, -E, g), s.lineTo(O.cx + O.x01, O.cy + O.y01), E < T ? s.arc(O.cx, O.cy, E, Z(O.y01, O.x01), Z(k.y01, k.x01), !g) : (s.arc(O.cx, O.cy, E, Z(O.y01, O.x01), Z(O.y11, O.x11), !g), s.arc(0, 0, d, Z(O.cy + O.y11, O.cx + O.x11), Z(k.cy + k.y11, k.cx + k.x11), g), s.arc(k.cx, k.cy, E, Z(k.y11, k.x11), Z(k.y01, k.x01), !g))) : s.arc(0, 0, d, b, y, g);
		}
		if (s.closePath(), l) return s = null, l + "" || null;
	}
	return l.centroid = function() {
		var n = (+e.apply(this, arguments) + +t.apply(this, arguments)) / 2, r = (+i.apply(this, arguments) + +a.apply(this, arguments)) / 2 - vs / 2;
		return [ms(r) * n, Q(r) * n];
	}, l.innerRadius = function(t) {
		return arguments.length ? (e = typeof t == "function" ? t : X(+t), l) : e;
	}, l.outerRadius = function(e) {
		return arguments.length ? (t = typeof e == "function" ? e : X(+e), l) : t;
	}, l.cornerRadius = function(e) {
		return arguments.length ? (n = typeof e == "function" ? e : X(+e), l) : n;
	}, l.padRadius = function(e) {
		return arguments.length ? (r = e == null ? null : typeof e == "function" ? e : X(+e), l) : r;
	}, l.startAngle = function(e) {
		return arguments.length ? (i = typeof e == "function" ? e : X(+e), l) : i;
	}, l.endAngle = function(e) {
		return arguments.length ? (a = typeof e == "function" ? e : X(+e), l) : a;
	}, l.padAngle = function(e) {
		return arguments.length ? (o = typeof e == "function" ? e : X(+e), l) : o;
	}, l.context = function(e) {
		return arguments.length ? (s = e ?? null, l) : s;
	}, l;
}
//#endregion
//#region node_modules/d3-shape/src/array.js
var Ms = Array.prototype.slice;
function Ns(e) {
	return typeof e == "object" && "length" in e ? e : Array.from(e);
}
//#endregion
//#region node_modules/d3-shape/src/point.js
function Ps(e) {
	return e[0];
}
function Fs(e) {
	return e[1];
}
//#endregion
//#region node_modules/d3-shape/src/curve/bump.js
var Is = class {
	constructor(e, t) {
		this._context = e, this._x = t;
	}
	areaStart() {
		this._line = 0;
	}
	areaEnd() {
		this._line = NaN;
	}
	lineStart() {
		this._point = 0;
	}
	lineEnd() {
		(this._line || this._line !== 0 && this._point === 1) && this._context.closePath(), this._line = 1 - this._line;
	}
	point(e, t) {
		switch (e = +e, t = +t, this._point) {
			case 0:
				this._point = 1, this._line ? this._context.lineTo(e, t) : this._context.moveTo(e, t);
				break;
			case 1: this._point = 2;
			default: this._x ? this._context.bezierCurveTo(this._x0 = (this._x0 + e) / 2, this._y0, this._x0, t, e, t) : this._context.bezierCurveTo(this._x0, this._y0 = (this._y0 + t) / 2, e, this._y0, e, t);
		}
		this._x0 = e, this._y0 = t;
	}
};
function Ls(e) {
	return new Is(e, !0);
}
//#endregion
//#region node_modules/d3-shape/src/link.js
function Rs(e) {
	return e.source;
}
function zs(e) {
	return e.target;
}
function Bs(e) {
	let t = Rs, n = zs, r = Ps, i = Fs, a = null, o = null, s = Cs(c);
	function c() {
		let c, l = Ms.call(arguments), u = t.apply(this, l), d = n.apply(this, l);
		if (a ?? (o = e(c = s())), o.lineStart(), l[0] = u, o.point(+r.apply(this, l), +i.apply(this, l)), l[0] = d, o.point(+r.apply(this, l), +i.apply(this, l)), o.lineEnd(), c) return o = null, c + "" || null;
	}
	return c.source = function(e) {
		return arguments.length ? (t = e, c) : t;
	}, c.target = function(e) {
		return arguments.length ? (n = e, c) : n;
	}, c.x = function(e) {
		return arguments.length ? (r = typeof e == "function" ? e : X(+e), c) : r;
	}, c.y = function(e) {
		return arguments.length ? (i = typeof e == "function" ? e : X(+e), c) : i;
	}, c.context = function(t) {
		return arguments.length ? (t == null ? a = o = null : o = e(a = t), c) : a;
	}, c;
}
function Vs() {
	return Bs(Ls);
}
//#endregion
//#region node_modules/d3-shape/src/offset/none.js
function Hs(e, t) {
	if ((o = e.length) > 1) for (var n = 1, r, i, a = e[t[0]], o, s = a.length; n < o; ++n) for (i = a, a = e[t[n]], r = 0; r < s; ++r) a[r][1] += a[r][0] = isNaN(i[r][1]) ? i[r][0] : i[r][1];
}
//#endregion
//#region node_modules/d3-shape/src/order/none.js
function Us(e) {
	for (var t = e.length, n = Array(t); --t >= 0;) n[t] = t;
	return n;
}
//#endregion
//#region node_modules/d3-shape/src/stack.js
function Ws(e, t) {
	return e[t];
}
function Gs(e) {
	let t = [];
	return t.key = e, t;
}
function Ks() {
	var e = X([]), t = Us, n = Hs, r = Ws;
	function i(i) {
		var a = Array.from(e.apply(this, arguments), Gs), o, s = a.length, c = -1, l;
		for (let e of i) for (o = 0, ++c; o < s; ++o) (a[o][c] = [0, +r(e, a[o].key, c, i)]).data = e;
		for (o = 0, l = Ns(t(a)); o < s; ++o) a[l[o]].index = o;
		return n(a, l), a;
	}
	return i.keys = function(t) {
		return arguments.length ? (e = typeof t == "function" ? t : X(Array.from(t)), i) : e;
	}, i.value = function(e) {
		return arguments.length ? (r = typeof e == "function" ? e : X(+e), i) : r;
	}, i.order = function(e) {
		return arguments.length ? (t = e == null ? Us : typeof e == "function" ? e : X(Array.from(e)), i) : t;
	}, i.offset = function(e) {
		return arguments.length ? (n = e ?? Hs, i) : n;
	}, i;
}
//#endregion
//#region node_modules/d3-zoom/src/constant.js
var qs = (e) => () => e;
//#endregion
//#region node_modules/d3-zoom/src/event.js
function Js(e, { sourceEvent: t, target: n, transform: r, dispatch: i }) {
	Object.defineProperties(this, {
		type: {
			value: e,
			enumerable: !0,
			configurable: !0
		},
		sourceEvent: {
			value: t,
			enumerable: !0,
			configurable: !0
		},
		target: {
			value: n,
			enumerable: !0,
			configurable: !0
		},
		transform: {
			value: r,
			enumerable: !0,
			configurable: !0
		},
		_: { value: i }
	});
}
//#endregion
//#region node_modules/d3-zoom/src/transform.js
function $(e, t, n) {
	this.k = e, this.x = t, this.y = n;
}
$.prototype = {
	constructor: $,
	scale: function(e) {
		return e === 1 ? this : new $(this.k * e, this.x, this.y);
	},
	translate: function(e, t) {
		return e === 0 & t === 0 ? this : new $(this.k, this.x + this.k * e, this.y + this.k * t);
	},
	apply: function(e) {
		return [e[0] * this.k + this.x, e[1] * this.k + this.y];
	},
	applyX: function(e) {
		return e * this.k + this.x;
	},
	applyY: function(e) {
		return e * this.k + this.y;
	},
	invert: function(e) {
		return [(e[0] - this.x) / this.k, (e[1] - this.y) / this.k];
	},
	invertX: function(e) {
		return (e - this.x) / this.k;
	},
	invertY: function(e) {
		return (e - this.y) / this.k;
	},
	rescaleX: function(e) {
		return e.copy().domain(e.range().map(this.invertX, this).map(e.invert, e));
	},
	rescaleY: function(e) {
		return e.copy().domain(e.range().map(this.invertY, this).map(e.invert, e));
	},
	toString: function() {
		return "translate(" + this.x + "," + this.y + ") scale(" + this.k + ")";
	}
};
var Ys = new $(1, 0, 0);
Xs.prototype = $.prototype;
function Xs(e) {
	for (; !e.__zoom;) if (!(e = e.parentNode)) return Ys;
	return e.__zoom;
}
//#endregion
//#region node_modules/d3-zoom/src/noevent.js
function Zs(e) {
	e.stopImmediatePropagation();
}
function Qs(e) {
	e.preventDefault(), e.stopImmediatePropagation();
}
//#endregion
//#region node_modules/d3-zoom/src/zoom.js
function $s(e) {
	return (!e.ctrlKey || e.type === "wheel") && !e.button;
}
function ec() {
	var e = this;
	return e instanceof SVGElement ? (e = e.ownerSVGElement || e, e.hasAttribute("viewBox") ? (e = e.viewBox.baseVal, [[e.x, e.y], [e.x + e.width, e.y + e.height]]) : [[0, 0], [e.width.baseVal.value, e.height.baseVal.value]]) : [[0, 0], [e.clientWidth, e.clientHeight]];
}
function tc() {
	return this.__zoom || Ys;
}
function nc(e) {
	return -e.deltaY * (e.deltaMode === 1 ? .05 : e.deltaMode ? 1 : .002) * (e.ctrlKey ? 10 : 1);
}
function rc() {
	return navigator.maxTouchPoints || "ontouchstart" in this;
}
function ic(e, t, n) {
	var r = e.invertX(t[0][0]) - n[0][0], i = e.invertX(t[1][0]) - n[1][0], a = e.invertY(t[0][1]) - n[0][1], o = e.invertY(t[1][1]) - n[1][1];
	return e.translate(i > r ? (r + i) / 2 : Math.min(0, r) || Math.max(0, i), o > a ? (a + o) / 2 : Math.min(0, a) || Math.max(0, o));
}
function ac() {
	var e = $s, t = ec, n = ic, r = nc, i = rc, a = [0, Infinity], o = [[-Infinity, -Infinity], [Infinity, Infinity]], s = 250, c = Br, l = ne("start", "zoom", "end"), u, d, f, p = 500, m = 150, h = 0, g = 10;
	function _(e) {
		e.property("__zoom", tc).on("wheel.zoom", w, { passive: !1 }).on("mousedown.zoom", T).on("dblclick.zoom", E).filter(i).on("touchstart.zoom", D).on("touchmove.zoom", O).on("touchend.zoom touchcancel.zoom", k).style("-webkit-tap-highlight-color", "rgba(0,0,0,0)");
	}
	_.transform = function(e, t, n, r) {
		var i = e.selection ? e.selection() : e;
		i.property("__zoom", tc), e === i ? i.interrupt().each(function() {
			S(this, arguments).event(r).start().zoom(null, typeof t == "function" ? t.apply(this, arguments) : t).end();
		}) : x(e, t, n, r);
	}, _.scaleBy = function(e, t, n, r) {
		_.scaleTo(e, function() {
			return this.__zoom.k * (typeof t == "function" ? t.apply(this, arguments) : t);
		}, n, r);
	}, _.scaleTo = function(e, r, i, a) {
		_.transform(e, function() {
			var e = t.apply(this, arguments), a = this.__zoom, s = i == null ? b(e) : typeof i == "function" ? i.apply(this, arguments) : i, c = a.invert(s), l = typeof r == "function" ? r.apply(this, arguments) : r;
			return n(y(v(a, l), s, c), e, o);
		}, i, a);
	}, _.translateBy = function(e, r, i, a) {
		_.transform(e, function() {
			return n(this.__zoom.translate(typeof r == "function" ? r.apply(this, arguments) : r, typeof i == "function" ? i.apply(this, arguments) : i), t.apply(this, arguments), o);
		}, null, a);
	}, _.translateTo = function(e, r, i, a, s) {
		_.transform(e, function() {
			var e = t.apply(this, arguments), s = this.__zoom, c = a == null ? b(e) : typeof a == "function" ? a.apply(this, arguments) : a;
			return n(Ys.translate(c[0], c[1]).scale(s.k).translate(typeof r == "function" ? -r.apply(this, arguments) : -r, typeof i == "function" ? -i.apply(this, arguments) : -i), e, o);
		}, a, s);
	};
	function v(e, t) {
		return t = Math.max(a[0], Math.min(a[1], t)), t === e.k ? e : new $(t, e.x, e.y);
	}
	function y(e, t, n) {
		var r = t[0] - n[0] * e.k, i = t[1] - n[1] * e.k;
		return r === e.x && i === e.y ? e : new $(e.k, r, i);
	}
	function b(e) {
		return [(+e[0][0] + +e[1][0]) / 2, (+e[0][1] + +e[1][1]) / 2];
	}
	function x(e, n, r, i) {
		e.on("start.zoom", function() {
			S(this, arguments).event(i).start();
		}).on("interrupt.zoom end.zoom", function() {
			S(this, arguments).event(i).end();
		}).tween("zoom", function() {
			var e = this, a = arguments, o = S(e, a).event(i), s = t.apply(e, a), l = r == null ? b(s) : typeof r == "function" ? r.apply(e, a) : r, u = Math.max(s[1][0] - s[0][0], s[1][1] - s[0][1]), d = e.__zoom, f = typeof n == "function" ? n.apply(e, a) : n, p = c(d.invert(l).concat(u / d.k), f.invert(l).concat(u / f.k));
			return function(e) {
				if (e === 1) e = f;
				else {
					var t = p(e), n = u / t[2];
					e = new $(n, l[0] - t[0] * n, l[1] - t[1] * n);
				}
				o.zoom(null, e);
			};
		});
	}
	function S(e, t, n) {
		return !n && e.__zooming || new C(e, t);
	}
	function C(e, n) {
		this.that = e, this.args = n, this.active = 0, this.sourceEvent = null, this.extent = t.apply(e, n), this.taps = 0;
	}
	C.prototype = {
		event: function(e) {
			return e && (this.sourceEvent = e), this;
		},
		start: function() {
			return ++this.active === 1 && (this.that.__zooming = this, this.emit("start")), this;
		},
		zoom: function(e, t) {
			return this.mouse && e !== "mouse" && (this.mouse[1] = t.invert(this.mouse[0])), this.touch0 && e !== "touch" && (this.touch0[1] = t.invert(this.touch0[0])), this.touch1 && e !== "touch" && (this.touch1[1] = t.invert(this.touch1[0])), this.that.__zoom = t, this.emit("zoom"), this;
		},
		end: function() {
			return --this.active === 0 && (delete this.that.__zooming, this.emit("end")), this;
		},
		emit: function(e) {
			var t = R(this.that).datum();
			l.call(e, this.that, new Js(e, {
				sourceEvent: this.sourceEvent,
				target: _,
				type: e,
				transform: this.that.__zoom,
				dispatch: l
			}), t);
		}
	};
	function w(t, ...i) {
		if (!e.apply(this, arguments)) return;
		var s = S(this, i).event(t), c = this.__zoom, l = Math.max(a[0], Math.min(a[1], c.k * 2 ** r.apply(this, arguments))), u = nn(t);
		if (s.wheel) (s.mouse[0][0] !== u[0] || s.mouse[0][1] !== u[1]) && (s.mouse[1] = c.invert(s.mouse[0] = u)), clearTimeout(s.wheel);
		else if (c.k === l) return;
		else s.mouse = [u, c.invert(u)], mi(this), s.start();
		Qs(t), s.wheel = setTimeout(d, m), s.zoom("mouse", n(y(v(c, l), s.mouse[0], s.mouse[1]), s.extent, o));
		function d() {
			s.wheel = null, s.end();
		}
	}
	function T(t, ...r) {
		if (f || !e.apply(this, arguments)) return;
		var i = t.currentTarget, a = S(this, r, !0).event(t), s = R(t.view).on("mousemove.zoom", d, !0).on("mouseup.zoom", p, !0), c = nn(t, i), l = t.clientX, u = t.clientY;
		on(t.view), Zs(t), a.mouse = [c, this.__zoom.invert(c)], mi(this), a.start();
		function d(e) {
			if (Qs(e), !a.moved) {
				var t = e.clientX - l, r = e.clientY - u;
				a.moved = t * t + r * r > h;
			}
			a.event(e).zoom("mouse", n(y(a.that.__zoom, a.mouse[0] = nn(e, i), a.mouse[1]), a.extent, o));
		}
		function p(e) {
			s.on("mousemove.zoom mouseup.zoom", null), sn(e.view, a.moved), Qs(e), a.event(e).end();
		}
	}
	function E(r, ...i) {
		if (e.apply(this, arguments)) {
			var a = this.__zoom, c = nn(r.changedTouches ? r.changedTouches[0] : r, this), l = a.invert(c), u = a.k * (r.shiftKey ? .5 : 2), d = n(y(v(a, u), c, l), t.apply(this, i), o);
			Qs(r), s > 0 ? R(this).transition().duration(s).call(x, d, c, r) : R(this).call(_.transform, d, c, r);
		}
	}
	function D(t, ...n) {
		if (e.apply(this, arguments)) {
			var r = t.touches, i = r.length, a = S(this, n, t.changedTouches.length === i).event(t), o, s, c, l;
			for (Zs(t), s = 0; s < i; ++s) c = r[s], l = nn(c, this), l = [
				l,
				this.__zoom.invert(l),
				c.identifier
			], a.touch0 ? !a.touch1 && a.touch0[2] !== l[2] && (a.touch1 = l, a.taps = 0) : (a.touch0 = l, o = !0, a.taps = 1 + !!u);
			u &&= clearTimeout(u), o && (a.taps < 2 && (d = l[0], u = setTimeout(function() {
				u = null;
			}, p)), mi(this), a.start());
		}
	}
	function O(e, ...t) {
		if (this.__zooming) {
			var r = S(this, t).event(e), i = e.changedTouches, a = i.length, s, c, l, u;
			for (Qs(e), s = 0; s < a; ++s) c = i[s], l = nn(c, this), r.touch0 && r.touch0[2] === c.identifier ? r.touch0[0] = l : r.touch1 && r.touch1[2] === c.identifier && (r.touch1[0] = l);
			if (c = r.that.__zoom, r.touch1) {
				var d = r.touch0[0], f = r.touch0[1], p = r.touch1[0], m = r.touch1[1], h = (h = p[0] - d[0]) * h + (h = p[1] - d[1]) * h, g = (g = m[0] - f[0]) * g + (g = m[1] - f[1]) * g;
				c = v(c, Math.sqrt(h / g)), l = [(d[0] + p[0]) / 2, (d[1] + p[1]) / 2], u = [(f[0] + m[0]) / 2, (f[1] + m[1]) / 2];
			} else if (r.touch0) l = r.touch0[0], u = r.touch0[1];
			else return;
			r.zoom("touch", n(y(c, l, u), r.extent, o));
		}
	}
	function k(e, ...t) {
		if (this.__zooming) {
			var n = S(this, t).event(e), r = e.changedTouches, i = r.length, a, o;
			for (Zs(e), f && clearTimeout(f), f = setTimeout(function() {
				f = null;
			}, p), a = 0; a < i; ++a) o = r[a], n.touch0 && n.touch0[2] === o.identifier ? delete n.touch0 : n.touch1 && n.touch1[2] === o.identifier && delete n.touch1;
			if (n.touch1 && !n.touch0 && (n.touch0 = n.touch1, delete n.touch1), n.touch0) n.touch0[1] = this.__zoom.invert(n.touch0[0]);
			else if (n.end(), n.taps === 2 && (o = nn(o, this), Math.hypot(d[0] - o[0], d[1] - o[1]) < g)) {
				var s = R(this).on("dblclick.zoom");
				s && s.apply(this, arguments);
			}
		}
	}
	return _.wheelDelta = function(e) {
		return arguments.length ? (r = typeof e == "function" ? e : qs(+e), _) : r;
	}, _.filter = function(t) {
		return arguments.length ? (e = typeof t == "function" ? t : qs(!!t), _) : e;
	}, _.touchable = function(e) {
		return arguments.length ? (i = typeof e == "function" ? e : qs(!!e), _) : i;
	}, _.extent = function(e) {
		return arguments.length ? (t = typeof e == "function" ? e : qs([[+e[0][0], +e[0][1]], [+e[1][0], +e[1][1]]]), _) : t;
	}, _.scaleExtent = function(e) {
		return arguments.length ? (a[0] = +e[0], a[1] = +e[1], _) : [a[0], a[1]];
	}, _.translateExtent = function(e) {
		return arguments.length ? (o[0][0] = +e[0][0], o[1][0] = +e[1][0], o[0][1] = +e[0][1], o[1][1] = +e[1][1], _) : [[o[0][0], o[0][1]], [o[1][0], o[1][1]]];
	}, _.constrain = function(e) {
		return arguments.length ? (n = e, _) : n;
	}, _.duration = function(e) {
		return arguments.length ? (s = +e, _) : s;
	}, _.interpolate = function(e) {
		return arguments.length ? (c = e, _) : c;
	}, _.on = function() {
		var e = l.on.apply(l, arguments);
		return e === l ? _ : e;
	}, _.clickDistance = function(e) {
		return arguments.length ? (h = (e = +e) * e, _) : Math.sqrt(h);
	}, _.tapDistance = function(e) {
		return arguments.length ? (g = +e, _) : g;
	}, _;
}
//#endregion
//#region src/Settings.ts
var oc = class {
	constructor() {
		this.width = 800, this.height = 800, this.enableTooltips = !0, this.tooltipContainer = null;
	}
}, sc;
(function(e) {
	function t(e) {
		return e < .5 ? 4 * e * e * e : 1 - (-2 * e + 2) ** 3 / 2;
	}
	e.easeInEaseOutCubic = t;
	function n(e) {
		return e * e * e;
	}
	e.easeInCubic = n;
	function r(e) {
		return 1 - (1 - e) ** 3;
	}
	e.easeOutCubic = r;
	function i(e) {
		let t = 2 * Math.PI / 4.5;
		return e === 0 ? 0 : e === 1 ? 1 : e < .5 ? -(2 ** (20 * e - 10) * Math.sin((20 * e - 11.125) * t)) / 2 : 2 ** (-20 * e + 10) * Math.sin((20 * e - 11.125) * t) / 2 + 1;
	}
	e.easeInEaseOutElastic = i;
	function a(e) {
		let t = 2 * Math.PI / 3;
		return e === 0 ? 0 : e === 1 ? 1 : -(2 ** (10 * e - 10)) * Math.sin((e * 10 - 10.75) * t);
	}
	e.easeInElastic = a;
	function o(e) {
		let t = 2 * Math.PI / 3;
		return e === 0 ? 0 : e === 1 ? 1 : 2 ** (-10 * e) * Math.sin((e * 10 - .75) * t) + 1;
	}
	e.easeOutElastic = o;
})(sc ||= {});
//#endregion
//#region src/visualizations/heatmap/cluster/TreeNode.ts
var cc = class e {
	static {
		this.currentID = 0;
	}
	constructor(t, n, r, i, a) {
		this._parent = t, this._leftChild = n, this._rightChild = r, this.values = i, this.height = a, this.id = e.currentID, e.currentID++;
	}
	get parent() {
		return this._parent;
	}
	set parent(e) {
		this._parent = e;
	}
	get leftChild() {
		return this._leftChild;
	}
	set leftChild(e) {
		this._leftChild = e;
	}
	get rightChild() {
		return this._rightChild;
	}
	set rightChild(e) {
		this._rightChild = e;
	}
	toNewick(e) {
		let t = "";
		return !this.leftChild && !this.rightChild ? e(this.values[0].id) + ":" + this.height : (t += "(", this.leftChild && (t += this.leftChild.toNewick(e) + ","), this.rightChild && (t += this.rightChild.toNewick(e)), t += ")" + this.id + ":" + this.height, t);
	}
	toGraphViz(e) {
		let t = "digraph dendrogram {\n", n = "", r = "", i = [this];
		for (; i.length > 0;) {
			let t = i.shift();
			if (!t) break;
			!t.leftChild && !t.rightChild ? n += `    ${t.id} [label="${e(t.values[0].id)}"];\n` : n += `    ${t.id} [label="${t.id}"];\n`, t.leftChild && (r += `    ${t.id} -> ${t.leftChild.id};\n`, i.push(t.leftChild)), t.rightChild && (r += `    ${t.id} -> ${t.rightChild.id};\n`, i.push(t.rightChild));
		}
		return t += n + r + "}", t;
	}
}, lc = class {
	constructor(e, t, n) {
		this.elements = e, this.index = t, this.treeNode = n;
	}
	merge(e, t) {
		this.elements.push(...e.elements);
		let n = new cc(null, this.treeNode, e.treeNode, this.elements.slice(), t);
		this.treeNode.parent = n, e.treeNode.parent = n, this.treeNode = n;
	}
}, uc = class {
	constructor(e) {
		this.metric = e;
	}
	cluster(e) {
		if (cc.currentID = 0, e.length < 1) return new cc(null, null, null, [], 0);
		let t = /* @__PURE__ */ new Map(), n = [];
		for (let r = 0; r < e.length; r++) {
			let i = e[r].values;
			t.set(r, new lc([e[r]], r, new cc(null, null, null, [e[r]], 0))), n.push(i);
		}
		let r = this.metric.getDistance(n), i = 0;
		for (; i != r.length - 1;) {
			let e = Infinity, n = -1, a = -1;
			for (let i of t.keys()) for (let o of t.keys()) i > o && r[i][o] < e && (e = r[i][o], n = i, a = o);
			let o = t.get(n), s = t.get(a), c = e / 2;
			if (!o || !s) throw "At least one cluster is invalid!";
			let l = this.copyDistanceMatrix(r);
			for (let e of t.keys()) if (e != n && e != a) {
				let t;
				t = e > n ? r[e][n] : r[n][e];
				let i;
				i = e > a ? r[e][a] : r[a][e];
				let c = (o.elements.length * t + s.elements.length * i) / (o.elements.length + s.elements.length);
				e > n ? l[e][n] = c : l[n][e] = c;
			}
			r = l, o.merge(s, c), t.delete(a), ++i;
		}
		return t.values().next().value.treeNode;
	}
	copyDistanceMatrix(e) {
		let t = [];
		for (let n = 0; n < e.length; n++) {
			let r = [], i = e[n];
			for (let e = 0; e < i.length; e++) r.push(i[e]);
			t.push(r);
		}
		return t;
	}
}, dc = class {
	getDistance(e) {
		let t = [];
		for (let n = 0; n < e.length; n++) {
			let r = [];
			for (let t = 0; t <= n; t++) r.push(this.calculateEuclideanDistance(e[n], e[t]));
			t.push(r);
		}
		return t;
	}
	calculateEuclideanDistance(e, t) {
		if (e.length != t.length) throw "Euclidean distance can only be calculated for 2 equally sized input arrays!";
		let n = 0;
		for (let r = 0; r < e.length; r++) n += (t[r] - e[r]) ** 2;
		return Math.sqrt(n);
	}
}, fc = class {
	constructor() {
		this.nodeMinMap = /* @__PURE__ */ new Map();
	}
	reorder(e) {
		return this.nodeMinMap.clear(), this.sortMinimum(e);
	}
	sortMinimum(e) {
		if (!e.leftChild || !e.rightChild) return e;
		let t = e.leftChild, n = e.rightChild, r = !t.leftChild && !t.rightChild, i = !n.leftChild && !n.rightChild;
		if (r && i) this.nodeMinMap.set(e, e.height);
		else if (!r && i) {
			let n = this.sortMinimum(t);
			e.leftChild = n;
			let r = this.nodeMinMap.get(n);
			if (r === void 0) throw "The recursive call to sort the left subtree did not yield a minimum value.";
			this.nodeMinMap.set(e, Math.min(e.height, r));
		} else if (r && !i) {
			let r = this.sortMinimum(n);
			e.leftChild = r, e.rightChild = t;
			let i = this.nodeMinMap.get(r);
			if (i === void 0) throw "The recursive call to sort the right subtree did not yield a minimum value.";
			this.nodeMinMap.set(e, Math.min(e.height, i));
		} else {
			let r = this.sortMinimum(t), i = this.sortMinimum(n), a = this.nodeMinMap.get(r), o = this.nodeMinMap.get(i);
			if (a === void 0 || o === void 0) throw "One of the recursive calls to sort a subtree did not yield a minimum value.";
			a <= o ? (e.leftChild = r, e.rightChild = i) : (e.leftChild = i, e.rightChild = r), this.nodeMinMap.set(e, Math.min(e.height, a, o));
		}
		return e;
	}
}, pc = class {
	constructor() {
		this.padding = {
			top: 8,
			right: 0,
			bottom: 0,
			left: 0
		}, this.title = "", this.titleFontSize = 14, this.width = 240, this.height = 12, this.labelFontSize = 12, this.ticks = 5, this.tickFormat = (e) => e.toFixed(2);
	}
}, mc = class extends oc {
	constructor(...e) {
		super(...e), this.initialTextWidth = 100, this.initialTextHeight = 100, this.squarePadding = 2, this.visualizationTextPadding = 4, this.fontSize = 14, this.labelColor = "#404040", this.highlightSelection = !0, this.highlightFontSize = 16, this.highlightFontColor = "black", this.className = "heatmap", this.animationsEnabled = !0, this.animationDuration = 2e3, this.transition = sc.easeInEaseOutCubic, this.minColor = "#EEEEEE", this.maxColor = "#1565C0", this.colorBuckets = 50, this.dendrogramEnabled = !1, this.dendrogramWidth = 100, this.dendrogramLineWidth = 1, this.dendrogramColor = "#404040", this.enableLegend = !1, this.legend = new pc(), this.clusteringAlgorithm = new uc(new dc()), this.reorderer = new fc(), this.getTooltip = (e, t, n) => `
            <style>
                .unipept-tooltip {
                    padding: 10px;
                    border-radius: 5px; 
                    background: rgba(0, 0, 0, 0.8); 
                    color: #fff;
                }
                
                .unipept-tooltip div, .unipept-tooltip a {
                    font-family: Roboto, 'Helvetica Neue', Helvetica, Arial, sans-serif;
                }
                
                .unipept-tooltip div {
                    font-weight: bold;
                }
            </style>
            <div class="unipept-tooltip">
                <div>
                    ${this.getTooltipTitle(e, t, n)}
                </div>
                <a>
                    ${this.getTooltipText(e)}
                </a>
            </div>
        `, this.getTooltipTitle = (e, t, n) => `${n.name ? n.name : ""}${n.name ? " and " : ""}${t.name ? t.name : ""}`, this.getTooltipText = (e) => `Similarity: ${(e.value * 100).toFixed(2)}%`;
	}
}, hc = class {
	constructor(e, t) {
		this.values = e, this.id = t;
	}
}, gc = class {
	preprocessFeatures(e) {
		return Object.entries(e).map(([e, t]) => ({
			name: t,
			idx: Number.parseInt(e)
		}));
	}
	preprocessValues(e, t, n, r) {
		let i = us().domain([0, 1]).range(this.computeColorPalette(t, n, r));
		return Object.entries(e).map(([e, t]) => Object.entries(t).map(([t, n]) => {
			if (typeof n == "number") {
				let r = i(n);
				if (r === void 0) throw Error("Invalid heatmap value given: " + n);
				return {
					value: n,
					rowId: Number.parseInt(e),
					columnId: Number.parseInt(t),
					color: r
				};
			}
			return n;
		}));
	}
	computeColorPalette(e, t, n) {
		let r = Vr(tr(e), tr(t));
		return Y().domain([0, 1]).range([0, 1]).ticks(n).map((e) => r(e));
	}
	orderPerColor(e) {
		let t = /* @__PURE__ */ new Map();
		for (let n = 0; n < e.length; n++) for (let r = 0; r < e[n].length; r++) {
			let i = e[n][r].color;
			t.has(i) || t.set(i, []), t.get(i)?.push([n, r]);
		}
		return t;
	}
}, _c = "tip", vc = "data-unipept-tooltip", yc = "data-unipept-tooltip-style", bc = class e {
	static {
		this.instances = /* @__PURE__ */ new WeakMap();
	}
	constructor(e) {
		this.container = e, this.element = null, this.visible = !1;
	}
	static create(t, n = null) {
		let r = e.resolveContainer(t, n), i = e.instances.get(r);
		return i || (i = new e(r), e.instances.set(r, i)), i;
	}
	static hideFor(t, n = null) {
		let r = e.findContainer(t, n);
		r && e.instances.get(r)?.hide();
	}
	show(e, t) {
		let n = this.ensureElement();
		n.innerHTML = t, this.moveTo(n, e), this.visible || (this.visible = !0, this.setVisibility(n, !0));
	}
	move(e) {
		this.moveTo(this.ensureElement(), e);
	}
	hide() {
		this.visible && this.element && (this.visible = !1, this.setVisibility(this.element, !1));
	}
	ensureElement() {
		return this.element && this.element.parentElement === this.container ? this.element : (this.element = this.container.querySelector(`:scope > .${_c}[${vc}]`) ?? this.createElement(), this.visible = !1, this.element);
	}
	moveTo(e, t) {
		e.style.top = `${t.clientY + 10}px`, e.style.left = `${t.clientX + 10}px`;
	}
	setVisibility(e, t) {
		e.style.visibility = t ? "visible" : "hidden", !(e.popover !== "manual" || !e.isConnected) && (t ? e.showPopover() : e.hidePopover());
	}
	createElement() {
		let t = this.container.ownerDocument.createElement("div");
		return t.className = _c, t.setAttribute(vc, ""), "popover" in t && (t.popover = "manual", e.installPopoverStyle(this.container.ownerDocument)), t.style.cssText = "position: fixed; inset: auto; z-index: 10; pointer-events: none; visibility: hidden;", this.container.appendChild(t), t;
	}
	static resolveContainer(t, n) {
		let r = e.findContainer(t, n);
		if (!r) throw Error(`No element matches the tooltip container selector "${n}".`);
		return r;
	}
	static findContainer(e, t) {
		return t ? typeof t == "string" ? e.ownerDocument.querySelector(t) : t : e.ownerDocument.body;
	}
	static installPopoverStyle(e) {
		if (e.querySelector(`style[${yc}]`)) return;
		let t = e.createElement("style");
		t.setAttribute(yc, ""), t.appendChild(e.createTextNode(`
            :where(.${_c}[popover]) {
                margin: 0;
                border: 0;
                padding: 0;
                background: transparent;
                color: inherit;
                overflow: visible;
            }
        `)), e.head.appendChild(t);
	}
}, xc = class {
	static clear(e, t = null) {
		bc.hideFor(e, t), e.innerHTML = "";
	}
}, Sc = class {
	constructor(e) {
		this.context = e;
	}
	renderLine(e, t, n, r, i, a) {
		this.context.lineWidth = i, this.context.moveTo(e, t), this.context.lineTo(n, r), this.context.strokeStyle = a, this.context.stroke();
	}
}, Cc = "'Helvetica Neue', Helvetica, Arial, sans-serif", wc = "rgba(0, 0, 0, 0.2)", Tc = 6, Ec = 8, Dc = 4, Oc = 1, kc = class {
	constructor(e, t, n, r, i = new mc()) {
		this.tooltip = null, this.legendElement = null, this.navigating = !1, this.highlightedRow = -1, this.highlightedColumn = -1, this.animatingRows = !1, this.animatingCols = !1, this.clusteredHorizontal = !1, this.clusteredVertical = !1, this.lastZoomStatus = {
			k: 1,
			x: 0,
			y: 0
		}, this.settings = this.fillOptions(i), this.element = e;
		let a = new gc();
		this.rows = a.preprocessFeatures(n), this.columns = a.preprocessFeatures(r), this.values = a.preprocessValues(t, this.settings.minColor, this.settings.maxColor, this.settings.colorBuckets), this.valuesPerColor = a.orderPerColor(this.values), this.settings.enableTooltips && (this.tooltip = bc.create(this.element, this.settings.tooltipContainer)), this.pixelRatio = window.devicePixelRatio || 1, this.originalViewPort = {
			xTop: 0,
			yTop: 0,
			xBottom: this.settings.width,
			yBottom: this.gridHeight
		}, this.currentViewPort = this.originalViewPort, this.textWidth = this.settings.initialTextWidth, this.textHeight = this.settings.initialTextHeight, xc.clear(this.element, this.settings.tooltipContainer), this.visElement = R(this.element).append("canvas").attr("width", this.pixelRatio * this.settings.width).attr("height", this.pixelRatio * this.gridHeight).attr("style", this.canvasStyle()).on("mouseover", (e) => this.tooltipMove(e)).on("mousemove", (e) => this.tooltipMove(e)).on("mouseout", (e) => this.tooltipMove(e)).on("click", (e) => this.click(e)), this.context = this.visElement.node().getContext("2d"), this.context.scale(this.pixelRatio, this.pixelRatio);
		let o = ac().extent([[0, 0], [this.settings.width, this.gridHeight]]).scaleExtent([.25, 12]).on("start", () => {
			this.navigating = !0, this.hideTooltip();
		}).on("zoom", (e) => {
			this.zoomed(e.transform);
		}).on("end", () => {
			this.navigating = !1;
		});
		this.visElement.call(o), this.settings.enableLegend && (this.legendElement = R(this.element).append("svg"), this.renderLegend()), this.computeClusterRoots(), this.redraw();
	}
	fillOptions(e = void 0) {
		let t = new mc();
		return Object.assign(t, e), t.legend = this.fillGroup(new pc(), e?.legend), t;
	}
	fillGroup(e, t) {
		let n = Object.assign({}, e.padding, t?.padding);
		return Object.assign(e, t, { padding: n });
	}
	reset() {
		this.redraw();
	}
	async cluster(e = "all") {
		let t = this.settings.animationsEnabled ? this.settings.animationDuration / 2 : 0, n = (e, n) => new Promise((r) => {
			let i, a = (o) => {
				i === void 0 && (i = o);
				let s = o - i, c = this.settings.transition(s / t);
				this.redraw(e, n, c), s < t ? requestAnimationFrame(a) : r();
			};
			requestAnimationFrame(a);
		}), r = new gc(), i = Array.from(Array(this.rows.length).keys()), a = Array(i.length);
		if ((e === "all" || e === "rows") && !this.clusteredVertical) {
			this.clusteredVertical = !0, i = this.determineOrder(this.rowClusterRoot);
			for (let [e, t] of Object.entries(i)) a[t] = Number.parseInt(e);
			let e = Array.from(Array(this.columns.length).keys());
			this.animatingRows = !0, await n(a, e), this.animatingRows = !1;
			let t = [];
			for (let e of i) t.push(this.values[e]);
			let o = [];
			for (let e of i) o.push(this.rows[e]);
			this.rows = o, this.values = t, this.valuesPerColor = r.orderPerColor(this.values);
		}
		let o = Array.from(Array(this.columns.length).keys()), s = Array(o.length);
		if ((e === "all" || e === "columns") && !this.clusteredHorizontal) {
			this.clusteredHorizontal = !0, o = this.determineOrder(this.colClusterRoot);
			for (let [e, t] of Object.entries(o)) s[t] = Number.parseInt(e);
			let e = Array.from(Array(this.rows.length).keys());
			this.animatingCols = !0, await n(e, s), this.animatingCols = !1;
			let t = [];
			for (let n of e) {
				let e = [];
				for (let t of o) e.push(this.values[n][t]);
				t.push(e);
			}
			let i = [];
			for (let e of o) i.push(this.columns[e]);
			this.columns = i, this.values = t, this.valuesPerColor = r.orderPerColor(this.values);
		}
		this.redraw();
	}
	computeClusterRoots() {
		let e = this.settings.clusteringAlgorithm, t = this.settings.reorderer, n = this.rows.map((e, t) => new hc(this.values[t].filter((t) => t.rowId == e.idx).map((e) => e.value), e.idx));
		this.rowClusterRoot = t.reorder(e.cluster(n)), this.verticalNodesPerDepth = this.bfsNodesPerDepth(this.rowClusterRoot);
		let r = this.columns.map((e, t) => new hc(this.values.map((e) => e[t].value), e.idx));
		this.colClusterRoot = t.reorder(e.cluster(r)), this.horizontalNodesPerDepth = this.bfsNodesPerDepth(this.colClusterRoot);
	}
	resize(e, t) {
		this.settings.width = e, this.settings.height = t, this.visElement.attr("height", this.pixelRatio * this.gridHeight), this.visElement.attr("width", this.pixelRatio * e), this.visElement.attr("style", this.canvasStyle()), this.context.scale(this.pixelRatio, this.pixelRatio), this.originalViewPort = {
			xTop: 0,
			yTop: 0,
			xBottom: e,
			yBottom: this.gridHeight
		}, this.renderLegend(), this.zoomed(this.lastZoomStatus);
	}
	get gridHeight() {
		return Math.max(Oc, this.settings.height - this.legendHeight());
	}
	canvasStyle() {
		return `${this.settings.enableLegend ? "display: block; " : ""}width: ${this.settings.width}px; height: ${this.gridHeight}px`;
	}
	legendHeight() {
		if (!this.settings.enableLegend) return 0;
		let e = this.settings.legend, t = e.title ? e.titleFontSize + Tc : 0;
		return e.padding.top + t + e.height + Ec + e.labelFontSize + e.padding.bottom;
	}
	legendBarWidth(e) {
		let t = this.settings.legend;
		return Math.max(1, Math.min(t.width, e - t.padding.left - t.padding.right));
	}
	computeLegendShapes(e) {
		let t = this.settings.legend, n = [], r = t.title ? t.titleFontSize + Tc : 0;
		t.title && n.push({
			type: "text",
			x: 0,
			y: 0,
			fontSize: t.titleFontSize,
			anchor: "start",
			content: t.title
		});
		let i = new gc().computeColorPalette(this.settings.minColor, this.settings.maxColor, this.settings.colorBuckets), a = e / i.length;
		for (let [o, s] of i.entries()) {
			let i = o * a;
			n.push({
				type: "rect",
				x: i,
				y: r,
				width: Math.min(a + 1, e - i),
				height: t.height,
				fill: s
			});
		}
		n.push({
			type: "rect",
			x: .5,
			y: r + .5,
			width: e - 1,
			height: t.height - 1,
			fill: "none",
			stroke: wc
		});
		let o = Y().domain([0, 1]).range([0, e]), s = o.ticks(t.ticks), c = r + t.height;
		for (let [e, r] of s.entries()) {
			let i = o(r);
			n.push({
				type: "line",
				x1: i,
				y1: c,
				x2: i,
				y2: c + Dc,
				stroke: this.settings.labelColor
			}), n.push({
				type: "text",
				x: i,
				y: c + Ec,
				fontSize: t.labelFontSize,
				anchor: e === 0 ? "start" : e === s.length - 1 ? "end" : "middle",
				content: t.tickFormat(r)
			});
		}
		return n;
	}
	renderLegend() {
		if (!this.legendElement) return;
		let e = this.settings.legend, t = this.legendHeight();
		this.legendElement.attr("width", this.settings.width).attr("height", t).attr("style", `display: block; width: ${this.settings.width}px; height: ${t}px`), this.legendElement.selectAll("*").remove();
		let n = this.legendElement.append("g").attr("class", "legend").attr("transform", `translate(${e.padding.left}, ${e.padding.top})`).attr("font-family", Cc).attr("fill", this.settings.labelColor);
		for (let e of this.computeLegendShapes(this.legendBarWidth(this.settings.width))) e.type === "rect" ? n.append("rect").attr("x", e.x).attr("y", e.y).attr("width", e.width).attr("height", e.height).attr("fill", e.fill).attr("stroke", e.stroke ?? null) : e.type === "line" ? n.append("line").attr("x1", e.x1).attr("y1", e.y1).attr("x2", e.x2).attr("y2", e.y2).attr("stroke", e.stroke) : n.append("text").attr("x", e.x).attr("y", e.y).attr("font-size", e.fontSize).attr("text-anchor", e.anchor).attr("dominant-baseline", "hanging").text(e.content);
	}
	legendShapesToSVG(e) {
		return e.map((e) => {
			if (e.type === "rect") {
				let t = e.stroke ? ` stroke="${e.stroke}"` : "";
				return `<rect x="${e.x}" y="${e.y}" width="${e.width}" height="${e.height}" fill="${e.fill}"${t}></rect>`;
			}
			return e.type === "line" ? `<line x1="${e.x1}" y1="${e.y1}" x2="${e.x2}" y2="${e.y2}" stroke="${e.stroke}"></line>` : `<text x="${e.x}" y="${e.y}" font-size="${e.fontSize}" text-anchor="${e.anchor}" dominant-baseline="hanging">${e.content}</text>`;
		}).join("\n");
	}
	toSVG(e = 14, t = 20, n = 2, r = 4) {
		let i = t, a = "";
		for (let [e, t] of this.valuesPerColor) for (let [r, o] of t) {
			let t = o * (i + n), s = r * (i + n);
			a += `
                    <rect width="${i}" height="${i}" fill="${e}" x="${t}" y="${s}"></rect>
                `;
		}
		let o = new OffscreenCanvas(1, 1).getContext("2d");
		o.font = `${e}px 'Helvetica Neue', Helvetica, Arial, sans-serif`;
		let s = i * this.columns.length + n * (this.columns.length - 1) + r, c = Math.max((i - e) / 2, 0), l = s;
		for (let t = 0; t < this.rows.length; t++) {
			let r = (i + n) * t + c;
			a += `
                <text 
                    x="${s}" 
                    y="${r}" 
                    font-size="${e}" 
                    dominant-baseline="hanging" 
                    fill="black"
                    font-family="'Helvetica Neue', Helvetica, Arial, sans-serif"
                >
                    ${this.rows[t].name}
                </text>
            `;
			let u = o.measureText(this.rows[t].name).width + s;
			u > l && (l = u);
		}
		let u = i * this.rows.length + n * (this.rows.length - 1) + r, d = u;
		for (let t = 0; t < this.columns.length; t++) {
			let r = (i + n) * t + c;
			a += `
                <text 
                    x="${r}" 
                    y="${u}" 
                    font-size="${e}" 
                    text-anchor="start" 
                    fill="black"
                    transform="rotate(90, ${r}, ${u})"
                    font-family="'Helvetica Neue', Helvetica, Arial, sans-serif"
                >
                    ${this.columns[t].name}
                </text>
            `;
			let s = o.measureText(this.columns[t].name).width + u;
			s > d && (d = s);
		}
		let f = "", p = 0;
		if (this.settings.enableLegend) {
			let e = this.settings.legend;
			p = this.legendHeight();
			let t = this.computeLegendShapes(this.legendBarWidth(l));
			f = `
                <g
                    class="legend"
                    transform="translate(${e.padding.left}, ${d + e.padding.top})"
                    font-family="${Cc}"
                    fill="${this.settings.labelColor}"
                >
                    ${this.legendShapesToSVG(t)}
                </g>
            `;
		}
		return `
            <svg xmlns="http://www.w3.org/2000/svg" width="${Math.ceil(l)}" height="${Math.ceil(d + p)}">
                ${a}
                ${f}
            </svg>
        `;
	}
	determineOrder(e) {
		return e.values.map((e) => e.id);
	}
	determineSquareWidth(e = this.currentViewPort, t = this.textWidth, n = this.textHeight) {
		let r = this.determineDendrogramWidth(), i = e.xBottom - e.xTop - r - this.columns.length * this.settings.squarePadding - t, a = e.yBottom - e.yTop - r - this.rows.length * this.settings.squarePadding - n, o = Math.max(1, i / this.columns.length), s = Math.max(1, a / this.rows.length);
		return Math.min(o, s);
	}
	determineDendrogramWidth() {
		return this.settings.dendrogramEnabled ? this.settings.dendrogramWidth * this.lastZoomStatus.k : 0;
	}
	computeTextStartX(e = this.currentViewPort, t = this.textWidth, n = this.textHeight) {
		return e.xTop + this.determineDendrogramWidth() + this.determineSquareWidth(e, t, n) * this.columns.length + this.settings.squarePadding * (this.columns.length - 1) + this.settings.visualizationTextPadding;
	}
	computeTextStartY(e = this.currentViewPort, t = this.textWidth, n = this.textHeight) {
		return e.yTop + this.determineDendrogramWidth() + this.determineSquareWidth(e, t, n) * this.rows.length + this.settings.squarePadding * (this.rows.length - 1) + this.settings.visualizationTextPadding;
	}
	zoomed({ k: e, x: t, y: n }) {
		this.lastZoomStatus = {
			k: e,
			x: t,
			y: n
		};
		let r = t + this.computeTextStartX(this.originalViewPort, this.settings.initialTextWidth, this.settings.initialTextHeight) * e, i = n + this.computeTextStartY(this.originalViewPort, this.settings.initialTextWidth, this.settings.initialTextHeight) * e, a = (t, n) => t > n ? n : e >= 1 ? Math.min(t, n) : Math.max(t, n);
		this.currentViewPort = {
			xTop: t + this.originalViewPort.xTop * e,
			yTop: n + this.originalViewPort.yTop * e,
			xBottom: a(t + this.originalViewPort.xBottom * e, this.originalViewPort.xBottom),
			yBottom: a(n + this.originalViewPort.yBottom * e, this.originalViewPort.yBottom)
		}, this.textWidth = this.currentViewPort.xBottom - r, this.textHeight = this.currentViewPort.yBottom - i, this.redraw();
	}
	redraw(e = Array.from(Array(this.rows.length).keys()), t = Array.from(Array(this.columns.length).keys()), n = -1) {
		this.redrawGrid(e, t, n), this.redrawRowTitles(e, n), this.redrawColumnTitles(t, n), this.redrawDendrogram(n);
	}
	redrawGrid(e, t, n) {
		n === -1 && (n = 0);
		let r = this.determineSquareWidth(), i = this.determineDendrogramWidth();
		this.context.clearRect(0, 0, this.settings.width, this.gridHeight);
		for (let [a, o] of this.valuesPerColor) {
			this.context.beginPath(), this.context.fillStyle = a;
			for (let [a, s] of o) {
				let o = this.currentViewPort.xTop + i + s * (r + this.settings.squarePadding), c = this.currentViewPort.yTop + i + a * (r + this.settings.squarePadding), l = this.currentViewPort.xTop + i + t[s] * (r + this.settings.squarePadding), u = this.currentViewPort.yTop + i + e[a] * (r + this.settings.squarePadding), d = l - o, f = u - c, p = o + d * n, m = c + f * n, h = p + (r + this.settings.squarePadding), g = m + (r + this.settings.squarePadding);
				h < 0 || p > this.settings.width || g < 0 || m > this.gridHeight || (this.settings.highlightSelection && a == this.highlightedRow && s == this.highlightedColumn && (this.context.save(), this.context.fillStyle = this.settings.maxColor, this.context.fillRect(p - this.settings.squarePadding, m - this.settings.squarePadding, r + 2 * this.settings.squarePadding, r + 2 * this.settings.squarePadding), this.context.restore()), this.context.fillRect(p, m, r, r));
			}
			this.context.closePath();
		}
	}
	ellipsizeString(e, t) {
		if (this.context.measureText(e).width > t) {
			let n = e.length, r = e.substr(0, n) + "...";
			for (; this.context.measureText(r).width > t && n > 0;) n--, r = e.substr(0, n) + "...";
			return n === 0 ? "" : r;
		}
		return e;
	}
	redrawRowTitles(e, t) {
		t === -1 && (t = 0);
		let n = this.determineSquareWidth(), r = this.determineDendrogramWidth(), i = Math.max(Math.floor((this.settings.fontSize + 12) / (n + this.settings.squarePadding)), 1), a = this.computeTextStartX(), o = Math.max((n - this.settings.fontSize) / 2, 0);
		this.context.save(), this.context.fillStyle = this.settings.labelColor, this.context.textBaseline = "top", this.context.textAlign = "start", this.context.font = `${this.settings.fontSize}px Arial, sans-serif`;
		for (let s = 0; s < this.rows.length; s += i) {
			let i = this.rows[s];
			this.settings.highlightSelection && s == this.highlightedRow && (this.context.save(), this.context.fillStyle = this.settings.highlightFontColor, this.context.font = `${this.settings.highlightFontSize}px 'Helvetica Neue', Helvetica, Arial, sans-serif`, o = Math.max((n - this.settings.highlightFontSize) / 2, 0));
			let c = this.currentViewPort.yTop + r + (n + this.settings.squarePadding) * s + o, l = c + (this.currentViewPort.yTop + r + (n + this.settings.squarePadding) * e[s] + o - c) * t;
			this.context.fillText(this.ellipsizeString(i.name, this.textWidth), a, l), this.settings.highlightSelection && s == this.highlightedRow && this.context.restore();
		}
		this.context.restore();
	}
	redrawColumnTitles(e, t) {
		t === -1 && (t = 0);
		let n = this.determineSquareWidth(), r = this.determineDendrogramWidth(), i = Math.max(Math.floor((this.settings.fontSize + 12) / (n + this.settings.squarePadding)), 1), a = this.computeTextStartY(), o = Math.max((n - this.settings.fontSize) / 2, 0);
		this.context.save(), this.context.rotate(90 * Math.PI / 180), this.context.fillStyle = this.settings.labelColor, this.context.textBaseline = "bottom", this.context.textAlign = "start", this.context.font = `${this.settings.fontSize}px Arial, sans-serif`;
		for (let s = 0; s < this.columns.length; s += i) {
			let i = this.columns[s];
			this.settings.highlightSelection && s == this.highlightedColumn && (this.context.save(), this.context.fillStyle = this.settings.highlightFontColor, this.context.font = `${this.settings.highlightFontSize}px 'Helvetica Neue', Helvetica, Arial, sans-serif`, o = Math.max((n - this.settings.highlightFontSize) / 2, 0));
			let c = -(this.currentViewPort.xTop + r + (n + this.settings.squarePadding) * s + o), l = c + (-(this.currentViewPort.xTop + r + (n + this.settings.squarePadding) * e[s] + o) - c) * t;
			this.context.fillText(this.ellipsizeString(i.name, this.textHeight), a, l), this.settings.highlightSelection && s == this.highlightedColumn && this.context.restore();
		}
		this.context.restore();
	}
	bfsNodesPerDepth(e) {
		let t = [], n = [];
		for (n.push([e, 0]); n.length > 0;) {
			let [e, r] = n.shift();
			t.length <= r && t.push([]), t[r].push(e), e.leftChild && n.push([e.leftChild, r + 1]), e.rightChild && n.push([e.rightChild, r + 1]);
		}
		return t;
	}
	redrawDendrogram(e) {
		this.settings.dendrogramEnabled && (this.redrawHorizontalDendrogram(e), this.redrawVerticalDendrogram(e));
	}
	computeDendrogramColor(e, t, n) {
		return n === -1 || !t ? e ? this.settings.dendrogramColor : "#d3d3d3" : Vr(tr("#d3d3d3"), tr(this.settings.dendrogramColor))(n);
	}
	redrawVerticalDendrogram(e) {
		this.context.save();
		let t = this.computeDendrogramColor(this.clusteredVertical, this.animatingRows, e), n = this.determineSquareWidth(), r = this.settings.dendrogramWidth * this.lastZoomStatus.k, i = new Sc(this.context), a = this.currentViewPort.yTop + r + n / 2, o = /* @__PURE__ */ new Map(), s = this.determineOrder(this.rowClusterRoot);
		for (let e = 0; e < s.length; e++) o.set(s[e], [this.currentViewPort.xTop + r, e * (n + this.settings.squarePadding) + a]);
		let c = r / this.rows.length, l = this.currentViewPort.xTop + r - c;
		for (let e = this.verticalNodesPerDepth.length - 1; e > 0; e--) for (let n = 0; n < this.verticalNodesPerDepth[e].length; n += 2) {
			let r = this.verticalNodesPerDepth[e][n], a = this.verticalNodesPerDepth[e][n + 1], s = r.parent, [u, d] = o.get(r.id), [f, p] = o.get(a.id);
			if (this.context.beginPath(), i.renderLine(u, d, l, d, this.settings.dendrogramLineWidth, t), i.renderLine(f, p, l, p, this.settings.dendrogramLineWidth, t), i.renderLine(l, d, l, p, this.settings.dendrogramLineWidth, t), this.context.closePath(), s) {
				let e = Math.min(d, p) + Math.abs(d - p) / 2;
				o.set(s.id, [l, e]);
			}
			l -= c;
		}
		if (!this.clusteredVertical) {
			this.context.rotate(-(90 * Math.PI) / 180), this.context.fillStyle = this.settings.labelColor;
			let e = 24 * this.lastZoomStatus.k;
			this.context.font = `${e}px 'Helvetica Neue', Helvetica, Arial, sans-serif`;
			let t = this.context.measureText("Click to cluster").width;
			this.context.fillText("Click to cluster", -(this.currentViewPort.yTop + r + this.rows.length * (n + this.settings.squarePadding) / 2) - t / 2, this.currentViewPort.xTop + r / 2 + e / 2);
		}
		this.context.restore();
	}
	redrawHorizontalDendrogram(e) {
		this.context.save();
		let t = this.computeDendrogramColor(this.clusteredHorizontal, this.animatingCols, e), n = this.determineSquareWidth(), r = this.settings.dendrogramWidth * this.lastZoomStatus.k, i = new Sc(this.context), a = this.currentViewPort.xTop + n / 2 + r, o = /* @__PURE__ */ new Map(), s = this.determineOrder(this.colClusterRoot);
		for (let e = 0; e < s.length; e++) o.set(s[e], [e * (n + this.settings.squarePadding) + a, this.currentViewPort.yTop + r]);
		let c = r / this.columns.length, l = this.currentViewPort.yTop + r - c;
		for (let e = this.horizontalNodesPerDepth.length - 1; e > 0; e--) for (let n = 0; n < this.horizontalNodesPerDepth[e].length; n += 2) {
			let r = this.horizontalNodesPerDepth[e][n], a = this.horizontalNodesPerDepth[e][n + 1], s = r.parent, [u, d] = o.get(r.id), [f, p] = o.get(a.id);
			if (this.context.beginPath(), i.renderLine(u, d, u, l, this.settings.dendrogramLineWidth, t), i.renderLine(f, p, f, l, this.settings.dendrogramLineWidth, t), i.renderLine(u, l, f, l, this.settings.dendrogramLineWidth, t), this.context.closePath(), s) {
				let e = Math.min(u, f) + Math.abs(u - f) / 2;
				o.set(s.id, [e, l]);
			}
			l -= c;
		}
		if (!this.clusteredHorizontal) {
			this.context.fillStyle = this.settings.labelColor;
			let e = 24 * this.lastZoomStatus.k;
			this.context.font = `${e}px 'Helvetica Neue', Helvetica, Arial, sans-serif`;
			let t = this.context.measureText("Click to cluster").width;
			this.context.fillText("Click to cluster", this.currentViewPort.xTop + r + this.columns.length * (n + this.settings.squarePadding) / 2 - t / 2, this.currentViewPort.yTop + r / 2 + e / 2);
		}
		this.context.restore();
	}
	findRowAndColForPosition(e, t) {
		let n = this.determineDendrogramWidth(), r = e - this.currentViewPort.xTop - n, i = t - this.currentViewPort.yTop - n, a = this.determineSquareWidth();
		return [Math.floor(i / (a + this.settings.squarePadding)), Math.floor(r / (a + this.settings.squarePadding))];
	}
	tooltipMove(e) {
		if (this.navigating) {
			this.hideTooltip();
			return;
		}
		let t = e.target.getBoundingClientRect(), [n, r] = this.findRowAndColForPosition(e.clientX - t.left, e.clientY - t.top);
		if (n < 0 || n >= this.rows.length || r < 0 || r >= this.columns.length) {
			this.hideTooltip(), this.highlightedRow = -1, this.highlightedColumn = -1, this.settings.highlightSelection && this.redraw();
			return;
		}
		this.highlightedRow = n, this.highlightedColumn = r, this.settings.highlightSelection && this.redraw(), this.settings.enableTooltips && this.tooltip && this.tooltip.show(e, this.settings.getTooltip(this.values[n][r], this.rows[n], this.columns[r]));
	}
	hideTooltip() {
		this.settings.enableTooltips && this.tooltip && this.tooltip.hide();
	}
	click(e) {
		if (!this.settings.dendrogramEnabled) return;
		let t = this.determineDendrogramWidth(), n = this.determineSquareWidth(), r = e.target.getBoundingClientRect(), i = e.clientX - r.left, a = e.clientY - r.top;
		if (i >= this.currentViewPort.xTop && i <= this.currentViewPort.xTop + t && a >= this.currentViewPort.yTop + t && a <= this.currentViewPort.yTop + t + this.rows.length * (n + this.settings.squarePadding)) {
			this.cluster("rows");
			return;
		}
		if (i >= this.currentViewPort.xTop + t && i <= this.currentViewPort.xTop + t + this.columns.length * (n + this.settings.squarePadding) && a >= this.currentViewPort.yTop && a <= this.currentViewPort.yTop + t) {
			this.cluster("columns");
			return;
		}
	}
}, Ac = class {
	getDistance(e) {
		let t = [];
		for (let n = 0; n < e.length; n++) {
			let r = [];
			for (let t = 0; t <= n; t++) r.push(this.getPearsonCorrelationBetween2Samples(e[n], e[t]));
			t.push(r);
		}
		return t;
	}
	getPearsonCorrelationBetween2Samples(e, t) {
		let n = (e, t) => e + t, r = e.reduce(n, 0) / e.length, i = t.reduce(n, 0) / t.length, a = 0, o = 0;
		for (let n = 0; n < e.length; n++) a += (e[n] - r) * (t[n] - i), o += Math.sqrt((e[n] - r) ** 2) * Math.sqrt((t[n] - i) ** 2);
		return 1 - a / o;
	}
}, jc = class {
	static {
		this.DEFAULT_COLORS = /* @__PURE__ */ "#f9f0ab.#e8e596.#f0e2a3.#ede487.#efd580.#f1cb82.#f1c298.#e8b598.#d5dda1.#c9d2b5.#aec1ad.#a7b8a8.#b49a3d.#b28647.#a97d32.#b68334.#d6a680.#dfad70.#a2765d.#9f6652.#b9763f.#bf6e5d.#af643c.#9b4c3f.#72659d.#8a6e9e.#8f5c85.#934b8b.#9d4e87.#92538c.#8b6397.#716084.#2e6093.#3a5988.#4a5072.#393e64.#aaa1cc.#e0b5c9.#e098b0.#ee82a2.#ef91ac.#eda994.#eeb798.#ecc099.#f6d5aa.#f0d48a.#efd95f.#eee469.#dbdc7f.#dfd961.#ebe378.#f5e351".split(".");
	}
	static {
		this.FIXED_COLORS = /* @__PURE__ */ "#1f77b4.#aec7e8.#ff7f0e.#ffbb78.#2ca02c.#98df8a.#d62728.#ff9896.#9467bd.#c5b0d5.#8c564b.#c49c94.#e377c2.#f7b6d2.#7f7f7f.#c7c7c7.#bcbd22.#dbdb8d.#17becf.#9edae5.#393b79.#5254a3.#6b6ecf.#9c9ede.#637939.#8ca252.#b5cf6b.#cedb9c.#8c6d31.#bd9e39.#e7ba52.#e7cb94.#843c39.#ad494a.#d6616b.#e7969c.#7b4173.#a55194.#ce6dbd.#de9ed6.#3182bd.#6baed6.#9ecae1.#c6dbef.#e6550d.#fd8d3c.#fdae6b.#fdd0a2.#31a354.#74c476.#a1d99b.#c7e9c0.#756bb1.#9e9ac8.#bcbddc.#dadaeb.#636363.#969696.#bdbdbd.#d9d9d9".split(".");
	}
	static {
		this.MATERIAL_DESIGN_COLORS = [
			"#ef5350",
			"#ec407a",
			"#ab47bc",
			"#7e57c2",
			"#5c6bc0",
			"#42a5f5",
			"#29b6f6",
			"#26c6da",
			"#26a69a",
			"#66bb6a",
			"#9ccc65",
			"#d4e157",
			"#ffee58",
			"#ffca28",
			"#ffa726",
			"#ff7043",
			"#8d6e63"
		];
	}
}, Mc = class {
	static stringHash(e) {
		return e.split("").reduce(function(e, t) {
			let n = (e << 5) - e + t.charCodeAt(0);
			return n & n;
		}, 0);
	}
}, Nc = class extends oc {
	constructor(...e) {
		super(...e), this.radius = 300, this.breadcrumbWidth = 200, this.className = "sunburst", this.useFixedColors = !1, this.colorPalette = jc.DEFAULT_COLORS, this.fixedColorPalette = jc.FIXED_COLORS, this.enableBreadcrumbs = !0, this.enableTooltips = !0, this.levels = 4, this.animationDuration = 1e3, this.rerootCallback = () => {}, this.fixedColorHash = (e) => Mc.stringHash(e.name), this.mouseIn = () => {}, this.mouseMove = () => {}, this.mouseOut = () => {}, this.getTooltip = (e) => `
            <style>
                .unipept-tooltip {
                    padding: 10px;
                    border-radius: 5px; 
                    background: rgba(0, 0, 0, 0.8); 
                    color: #fff;
                }
                
                .unipept-tooltip div, .unipept-tooltip a {
                    font-family: Roboto, 'Helvetica Neue', Helvetica, Arial, sans-serif;
                }
                
                .unipept-tooltip div {
                    font-weight: bold;
                }
            </style>
            <div class="unipept-tooltip">
                <div>
                    ${this.getTooltipTitle(e)}
                </div>
                <a>
                    ${this.getTooltipText(e)}
                </a>
            </div>
        `, this.getTooltipTitle = (e) => e.name, this.getTooltipText = (e) => `${e.count} hits`, this.getLabel = (e) => e.name === "empty" ? "" : e.name, this.getTitleText = this.getLabel;
	}
}, Pc = class t {
	static {
		this.idCounter = 0;
	}
	preprocessData(n) {
		let r = [];
		if (n.children) for (let e of n.children) r.push(this.preprocessData(e));
		return r.length > 0 && n.count !== 0 && r.push(new e(-1, "empty", [], n.count, n.selfCount)), new e(n.id || ++t.idCounter, n.name || "", r, n.count, n.selfCount, n.extra);
	}
}, Fc = class {
	static isParentOf(e, t, n) {
		if (t.depth >= n) return !1;
		let r = t;
		for (; r;) {
			if (r === e) return !0;
			r = r.parent;
		}
		return !1;
	}
}, Ic = class {
	static applyStyle(e, t, n) {
		e.classList.add(t);
		let r = e.ownerDocument, i = `style[data-unipept-style="${t}"]`, a = r.head.querySelector(i);
		a || (a = r.createElement("style"), a.setAttribute("data-unipept-style", t), r.head.appendChild(a)), a.textContent = n;
	}
}, Lc = class e {
	static getReadableColorFor(t) {
		return e.brightness(jn(t)) < 125 ? "#eee" : "#000";
	}
	static brightness({ r: e, g: t, b: n }) {
		return e * .299 + t * .587 + n * .114;
	}
}, Rc = class {
	constructor(e, t, n = new Nc()) {
		this.element = e, this.colorCounter = -1, this.blendedColors = /* @__PURE__ */ new Map(), this.currentMaxLevel = 4, this.arcData = [], this.textData = [], this.textKeys = /* @__PURE__ */ new Map(), this.renderGeneration = 0, this.previousRoot = null, this.previousMaxLevel = this.currentMaxLevel, this.settings = this.fillOptions(n);
		let r = new Pc().preprocessData(t);
		this.settings.enableTooltips && (this.tooltip = bc.create(this.element, this.settings.tooltipContainer)), this.currentMaxLevel = this.settings.levels, this.xScale = Y().range([0, 2 * Math.PI]), this.yScale = Y().domain([0, 1]).range([0, this.settings.radius]);
		let i = bo(r);
		i.sum((e) => e.children.length > 0 ? 0 : e.selfCount);
		let a = Mo();
		this.data = a(i).descendants(), this.data.forEach((e, t) => this.textKeys.set(e, t.toString())), this.arc = js().startAngle((e) => Math.max(0, Math.min(Math.PI * 2, this.xScale(e.x0)))).endAngle((e) => Math.max(0, Math.min(Math.PI * 2, this.xScale(e.x1)))).innerRadius((e) => Math.max(0, e.y0 ? this.yScale(e.y0) : e.y0)).outerRadius((e) => Math.max(0, this.yScale(e.y1) + 1)), this.initCss(), xc.clear(this.element, this.settings.tooltipContainer), this.breadCrumbs = R(this.element).append("div").attr("id", Math.floor(Math.random() * 2 ** 16) + "-breadcrumbs").attr("class", "sunburst-breadcrumbs").append("ul");
		let o = R(this.element).append("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg").attr("viewBox", `0 0 ${this.settings.width} ${this.settings.height}`).attr("width", this.settings.width).attr("height", this.settings.height).attr("overflow", "hidden").style("font-family", "'Helvetica Neue', Helvetica, Arial, sans-serif");
		o.append("style").attr("type", "text/css").html(".hidden{ visibility: hidden;}"), this.visGElement = o.append("g").attr("transform", "translate(" + this.settings.radius + "," + this.settings.radius + ")"), this.reset();
	}
	reset() {
		this.click(this.data[0]);
	}
	reroot(e, t = !0) {
		let n = this.data.find((t) => t.data.id === e);
		n && this.click(n, t);
	}
	fillOptions(e = void 0) {
		let t = new Nc();
		return Object.assign(t, e);
	}
	maxY(e) {
		return e.children ? Math.max(...e.children.map((e) => this.maxY(e))) : e.y1;
	}
	isBlended(e) {
		return e.name !== "empty" && !this.settings.useFixedColors && e.children.length > 0;
	}
	unblendedColor(e) {
		return e.name === "empty" ? "white" : this.settings.useFixedColors ? this.settings.fixedColorPalette[Math.abs(this.settings.fixedColorHash(e)) % this.settings.fixedColorPalette.length] : (e.extra.color || (e.extra.color = this.getColor()), e.extra.color);
	}
	blendedColor(e) {
		let t = this.blendedColors.get(e);
		if (t !== void 0) return t;
		let [n, r] = e.children.map((e) => this.isBlended(e) ? this.blendedColor(e) : Bn(this.unblendedColor(e))), i = e.children.length === 1 || e.children[1].name === "empty" ? Bn(n.h, n.s, n.l * .98) : Bn((n.h + r.h) / 2, (n.s + r.s) / 2, (n.l + r.l) / 2);
		return this.blendedColors.set(e, i), i;
	}
	color(e) {
		return this.isBlended(e) ? this.blendedColor(e).toString() : this.unblendedColor(e);
	}
	getColor() {
		return this.colorCounter = (this.colorCounter + 1) % this.settings.colorPalette.length, this.settings.colorPalette[this.colorCounter];
	}
	initCss() {
		let e = this.settings.className;
		Ic.applyStyle(this.element, e, `
.${e} {
    font-family: Roboto,'Helvetica Neue',Helvetica,Arial,sans-serif;
    width: ${this.settings.width + this.settings.breadcrumbWidth}px;
}
.${e} .sunburst-breadcrumbs {
    width: 176px;
    float: right;
    margin-right: 15px;
    margin-top: 10px;
    padding-left: 5px;
}
.${e} .sunburst-breadcrumbs ul {
    padding-left: 0;
    list-style: none;
}
.${e} .sunburst-breadcrumbs .crumb {
    margin-bottom: 5px;
    cursor: pointer;
}
.${e} .sunburst-breadcrumbs .crumb svg {
    float: left;
    margin-right: 3px;
}
.${e} .sunburst-breadcrumbs .crumb p {
    white-space: nowrap;
    text-overflow: ellipsis;
    overflow: hidden;
    margin: 0;
    font-size: 14px;
}
.${e} .sunburst-breadcrumbs .crumb .percentage {
    font-size: 11px;
}`);
	}
	arcTween(e, t) {
		let n = Math.min(this.maxY(e), e.y0 + t.settings.levels * (e.y1 - e.y0)), r = Tr(t.xScale.domain(), [e.x0, e.x1]), i = Tr(t.yScale.domain(), [e.y0, n]), a = Tr(t.yScale.range(), [e.y0 ? 20 : 0, t.settings.radius]);
		return (e) => (n) => (t.xScale.domain(r(n)), t.yScale.domain(i(n)).range(a(n)), t.arc(e));
	}
	tooltipIn(e, t) {
		this.settings.enableTooltips && this.tooltip && t.depth < this.currentMaxLevel && t.data.name !== "empty" && this.tooltip.show(e, this.settings.getTooltip(t.data));
	}
	tooltipMove(e, t) {
		this.settings.enableTooltips && this.tooltip && this.tooltip.move(e);
	}
	tooltipOut(e, t) {
		this.settings.enableTooltips && this.tooltip && this.tooltip.hide();
	}
	computeAvailableSpace(e) {
		return 2 * Math.max(0, this.yScale(e.y1) + 1) * Math.PI * (Math.max(0, Math.min(Math.PI * 2, this.xScale(e.x1)) - Math.max(0, Math.min(Math.PI * 2, this.xScale(e.x0)))) / (2 * Math.PI));
	}
	click(e, t = !0) {
		e.data.name === "empty" || this.previousRoot && this.previousRoot.data.id === e.data.id || (this.previousRoot = e, this.settings.enableBreadcrumbs && this.setBreadcrumbs(e), this.settings.rerootCallback && t && this.settings.rerootCallback(e.data), this.currentMaxLevel = e.depth + this.settings.levels, this.renderGeneration++, this.renderArcs(e), this.renderText(e));
	}
	async renderArcs(e) {
		let t = this.data.filter((t) => Fc.isParentOf(e, t, this.currentMaxLevel + 2));
		e.parent && t.push(e.parent);
		let n = t.filter((e) => !this.arcData.includes(e)), r = this.arcData.concat(...n);
		this.visGElement.selectAll("path").data([]).exit().remove(), this.path = this.visGElement.selectAll("path").data(r).enter().insert("path").attr("class", "arc").attr("id", (e, t) => "path-" + t).attr("d", this.arc).attr("fill-rule", "evenodd").style("fill", (e) => this.color(e.data)).attr("fill-opacity", (e) => e.depth >= this.previousMaxLevel ? .2 : 1).on("click", (e, t) => {
			t.depth < this.currentMaxLevel && this.click(t);
		}).on("mouseover", (e, t) => {
			this.settings.mouseIn(t.data, {
				x: e.clientX,
				y: e.clientY
			}), this.tooltipIn(e, t);
		}).on("mousemove", (e, t) => {
			this.settings.mouseMove(t.data, {
				x: e.clientX,
				y: e.clientY
			}), this.tooltipMove(e, t);
		}).on("mouseout", (e, t) => {
			this.settings.mouseOut(t.data), this.tooltipOut(e, t);
		}), await new Promise((t) => {
			this.path.transition().duration(this.settings.animationDuration).attrTween("d", this.arcTween(e, this)).attr("class", (e) => e.depth >= this.currentMaxLevel ? "arc toHide" : "arc").attr("fill-opacity", (e) => e.depth >= this.currentMaxLevel ? .2 : 1).on("end", () => {
				t();
			});
		}), this.previousMaxLevel = this.currentMaxLevel, this.arcData = t;
	}
	async renderText(e) {
		let t = this.renderGeneration, n = this.data.filter((t) => Fc.isParentOf(e, t, this.currentMaxLevel)), r = n.filter((e) => !this.textData.includes(e)), i = this.textData.concat(...r), a = e.parent ? i.indexOf(e.parent) : -1;
		a !== -1 && i.splice(a, 1);
		let o = this, s = typeof OffscreenCanvas < "u", c;
		s && (c = new OffscreenCanvas(1, 1).getContext("2d"), c.font = c.font = "16px 'Helvetica Neue', Helvetica, Arial, sans-serif");
		let l = this.visGElement.selectAll("text").data(i, (e) => this.textKeys.get(e));
		l.exit().remove(), this.text = l.enter().append("text").style("fill-opacity", 0).style("font-family", "font-family: Helvetica, 'Super Sans', sans-serif").style("pointer-events", "none").attr("dy", ".2em").merge(l).style("fill", (e) => Lc.getReadableColorFor(this.color(e.data))).style("visibility", null).text((e) => this.settings.getLabel(e.data)).style("font-size", function(e) {
			let t = s ? c.measureText(this.textContent).width : this.getComputedTextLength();
			return Math.floor(Math.min(o.settings.radius / o.settings.levels / t * 10 + 1, 12)) + "px";
		}).raise(), await new Promise((t) => {
			if (this.text.empty()) {
				t();
				return;
			}
			this.text.transition().duration(this.settings.animationDuration).attrTween("text-anchor", (e) => (t) => this.xScale(e.x0 + (e.x1 - e.x0) / 2) > Math.PI ? "end" : "start").attrTween("dx", (e) => (t) => this.xScale(e.x0 + (e.x1 - e.x0) / 2) > Math.PI ? "-4px" : "4px").attrTween("transform", (e) => (t) => {
				let n = this.xScale(e.x0 + (e.x1 - e.x0) / 2) * 180 / Math.PI - 90;
				return `rotate(${n})translate(${this.yScale(e.y0)})rotate(${n > 90 ? -180 : 0})`;
			}).styleTween("fill-opacity", function(e) {
				let t = Number.parseInt(R(this).style("font-size").replace("px", "")), n = Number.parseFloat(R(this).style("fill-opacity")) || 0;
				return (r) => o.computeAvailableSpace(e) > t ? (n + (1 - n) * r).toString() : "0";
			}).on("interrupt.resolve cancel.resolve", () => t()).on("end", function(n) {
				let r = o.computeAvailableSpace(n), i = R(this);
				i.style("visibility", r > Number.parseInt(i.style("font-size").replace("px", "")) && Fc.isParentOf(e, n, o.currentMaxLevel) ? "visible" : "hidden"), t();
			});
		}), t === this.renderGeneration && (this.textData = n);
	}
	setBreadcrumbs(e) {
		let t = [], n = e;
		for (; n;) t.push(n), n = n.parent;
		t.reverse().shift();
		let r = js().innerRadius(0).outerRadius(15).startAngle(0).endAngle((e) => 2 * Math.PI * e.data.count / e.parent.data.count);
		this.breadCrumbs.selectAll(".crumb").data(t).enter().append("li").on("click", (e, t) => {
			this.click(t.parent);
		}).attr("class", "crumb").style("opacity", "0").attr("title", (e) => this.settings.getTitleText(e.data)).html((e) => `
<p class='name'>${e.data.name}</p>
<p class='percentage'>${Math.round(100 * e.data.count / e.parent.data.count)}% of ${e.parent?.data.name}</p>`).insert("svg", ":first-child").attr("width", 30).attr("height", 30).append("path").attr("d", r).attr("transform", "translate(15, 15)").attr("fill", (e) => this.color(e.data)), this.breadCrumbs.selectAll(".crumb").transition().duration(this.settings.animationDuration).style("opacity", "1"), this.breadCrumbs.selectAll(".crumb").data(t).exit().transition().duration(this.settings.animationDuration).style("opacity", "0").remove();
	}
}, zc = class extends oc {
	constructor(...e) {
		super(...e), this.className = "treemap", this.levels = void 0, this.labelHeight = 10, this.colorRoot = "#104B7D", this.colorLeaf = "#fdffcc", this.colorBreadcrumbs = "#FF8F00", this.rerootCallback = () => {}, this.getBreadcrumbTooltip = (e) => e.name, this.getTooltip = (e) => `
            <style>
                .unipept-tooltip {
                    padding: 10px;
                    border-radius: 5px; 
                    background: rgba(0, 0, 0, 0.8); 
                    color: #fff;
                }
                
                .unipept-tooltip div, .unipept-tooltip a {
                    font-family: Roboto, 'Helvetica Neue', Helvetica, Arial, sans-serif;
                }
                
                .unipept-tooltip div {
                    font-weight: bold;
                }
            </style>
            <div class="unipept-tooltip">
                <div>
                    ${this.getTooltipTitle(e)}
                </div>
                <a>
                    ${this.getTooltipText(e)}
                </a>
            </div>
        `, this.getTooltipTitle = (e) => e.name, this.getTooltipText = (e) => `${e.count} hits`, this.getLabel = (e) => e.name, this.getLevel = (e) => e.depth;
	}
}, Bc = class t {
	static {
		this.idCounter = 0;
	}
	preprocessData(n) {
		let r = [];
		if (n.children) for (let e of n.children) r.push(this.preprocessData(e));
		return new e(n.id || ++t.idCounter, n.name || "", r, n.count, n.selfCount, n.extra);
	}
}, Vc = class {
	constructor(e, t, n = new zc()) {
		this.element = e, this.childParentRelations = /* @__PURE__ */ new Map(), this.nodeId = 0, this.settings = this.fillOptions(n), this.settings.enableTooltips && (this.tooltip = bc.create(this.element, this.settings.tooltipContainer)), this.initCss();
		let r = bo(new Bc().preprocessData(t));
		r.sum((e) => e.children.length > 0 ? 0 : e.count), r.sort((e, t) => t.value - e.value), this.partition = Ko(), this.partition.size([this.settings.width + 1, this.settings.height + 1]).paddingTop(this.settings.labelHeight), this.data = this.partition(r).descendants(), this.settings.levels || (this.settings.levels = this.data[0].height);
		for (let e of this.data) this.childParentRelations.set(e.data, e.parent?.data);
		this.currentRoot = this.data[0], this.colorScale = Y().domain([0, this.settings.levels]).range([this.settings.colorRoot, this.settings.colorLeaf]).interpolate(Vr), xc.clear(this.element, this.settings.tooltipContainer), this.breadCrumbs = R(this.element).append("div").attr("class", "breadcrumbs").style("position", "relative").style("width", this.settings.width + "px").style("height", "20px").style("background-color", this.settings.colorBreadcrumbs), this.treemap = R(this.element).append("div").style("position", "relative").style("width", this.settings.width + "px").style("height", this.settings.height + "px"), this.render(this.currentRoot);
	}
	resize(e, t) {
		this.settings.width = e, this.settings.height = t, this.partition.size([e + 1, t + 1]), this.breadCrumbs.style("width", this.settings.width + "px"), this.treemap.style("width", this.settings.width + "px"), this.treemap.style("height", this.settings.height + "px"), this.render(this.currentRoot, !1);
	}
	reroot(e, t = !0) {
		let n = this.data.find((t) => t.data.id === e);
		n && this.render(n, t);
	}
	reset() {
		this.render(this.data[0], !1);
	}
	fillOptions(e = void 0) {
		let t = new zc();
		return Object.assign(t, e);
	}
	initCss() {
		let e = this.settings.className;
		Ic.applyStyle(this.element, e, `
            .${e} {
                font-family: Arial,sans-serif;
            }
            .${e} .node {
                font-size: 9px;
                line-height: 10px;
                overflow: hidden;
                position: absolute;
                text-indent: 2px;
                text-align: center;
                text-overflow: ellipsis;
                cursor: pointer;
            }
            .${e} .node:hover {
                outline: 1px solid white;
            }
            .${e} .breadcrumbs {
                font-size: 11px;
                line-height: 20px;
                padding-left: 5px;
                font-weight: bold;
                color: white;
                box-sizing: border-box;
            }
            .full-screen .${e} .breadcrumbs {
                width: 100% !important;
            }
            .${e} .crumb {
                cursor: pointer;
            }
            .${e} .crumb .link:hover {
                text-decoration: underline;
            }
            .${e} .breadcrumbs .crumb + .crumb::before {
                content: " > ";
                cursor: default;
            }
        `);
	}
	render(e, t = !0) {
		this.currentRoot = e, this.setBreadcrumbs();
		let n = bo(e.data);
		n.sum((e) => e.children.length > 0 ? 0 : e.count), n.sort((e, t) => t.value - e.value);
		let r = this.treemap.selectAll(".node").data(this.partition(n).descendants(), (e) => e.data.id || (e.data.id = ++this.nodeId));
		r.enter().append("div").attr("class", "node").style("background", (e) => this.colorScale(this.settings.getLevel(e))).style("color", (e) => Lc.getReadableColorFor(this.colorScale(this.settings.getLevel(e)).toString())).style("left", "0px").style("top", "0px").style("width", "0px").style("height", "0px").text((e) => this.settings.getLabel(e.data)).on("click", (e, t) => this.render(t)).on("contextmenu", (e, t) => {
			e.preventDefault(), this.currentRoot.parent && this.render(this.currentRoot.parent);
		}).on("mouseover", (e, t) => this.tooltipIn(e, t)).on("mousemove", (e, t) => this.tooltipMove(e, t)).on("mouseout", (e, t) => this.tooltipOut(e, t)).merge(r).order().transition().call((e) => {
			e.style("left", (e) => e.x0 + "px"), e.style("top", (e) => e.y0 + "px"), e.style("width", (e) => Math.max(0, e.x1 - e.x0 - 1) + "px"), e.style("height", (e) => Math.max(0, e.y1 - e.y0 - 1) + "px");
		}), r.exit().remove(), t && this.settings.rerootCallback(this.currentRoot.data);
	}
	setBreadcrumbs() {
		let e = [], t = this.currentRoot.data;
		for (; t;) e.push(t), t = this.childParentRelations.get(t);
		e.reverse(), this.breadCrumbs.html(""), this.breadCrumbs.selectAll(".crumb").data(e).enter().append("span").attr("class", "crumb").attr("title", (e) => this.settings.getBreadcrumbTooltip(e)).html((e) => `<span class='link'>${e.name}</span>`).on("click", (e, t) => {
			this.render(this.data.filter((e) => e.data.id === t.id)[0]);
		});
	}
	tooltipIn(e, t) {
		this.settings.enableTooltips && this.tooltip && this.tooltip.show(e, this.settings.getTooltip(t.data));
	}
	tooltipMove(e, t) {
		this.settings.enableTooltips && this.tooltip && this.tooltip.move(e);
	}
	tooltipOut(e, t) {
		this.settings.enableTooltips && this.tooltip && this.tooltip.hide();
	}
}, Hc = Yo(fs), Uc = class extends oc {
	constructor(...e) {
		super(...e), this.minNodeSize = 2, this.maxNodeSize = 105, this.enableExpandOnClick = !0, this.enableAutoExpand = !1, this.autoExpandValue = .8, this.levelsToExpand = 2, this.enableRightClick = !0, this.enableInnerArcs = !0, this.enableLabels = !0, this.nodeDistance = 180, this.animationDuration = 500, this.colorProviderLevels = 1, this.nodeFillColor = (e) => e.isSelected() ? e.children.length > 0 ? e.getColor() || "#aaa" : "#fff" : "#aaa", this.nodeStrokeColor = (e) => e.isSelected() && e.getColor() || "#aaa", this.linkStrokeColor = (e) => e.source.data.isSelected() ? e.target.data.getColor() : "#aaa", this.colorProvider = (e) => Hc(e.name), this.getLabel = (e) => e.name, this.getTooltip = (e) => `
            <style>
                .unipept-tooltip {
                    padding: 10px;
                    border-radius: 5px; 
                    background: rgba(0, 0, 0, 0.8); 
                    color: #fff;
                }
                
                .unipept-tooltip div, .unipept-tooltip a {
                    font-family: Roboto, 'Helvetica Neue', Helvetica, Arial, sans-serif;
                }
                
                .unipept-tooltip div {
                    font-weight: bold;
                }
            </style>
            <div class="unipept-tooltip">
                <div>
                    ${this.getTooltipTitle(e)}
                </div>
                <a>
                    ${this.getTooltipText(e)}
                </a>
            </div>
        `, this.getTooltipTitle = (e) => e.name, this.getTooltipText = (e) => `${e.count} hits`;
	}
}, Wc = class {
	constructor(e = [], t) {
		this.data = e, this.comparator = t, this.heapify();
	}
	add(e) {
		this.data.push(e), this.bubbleUp(this.data.length - 1);
	}
	peek() {
		return this.data[0];
	}
	remove() {
		let e = this.data[0];
		return this.data.length > 1 ? (this.data[0] = this.data.pop(), this.sink(0)) : this.data.pop(), e;
	}
	clear() {
		this.data.splice(0, this.data.length);
	}
	size() {
		return this.data.length;
	}
	heapify() {
		for (let e = Math.floor((this.data.length - 2) / 2); e >= 0; e--) this.sink(e);
	}
	bubbleUp(e) {
		let t = this.data[e];
		for (; e > 0;) {
			let n = Math.floor((e - 1) / 2), r = this.data[n];
			if (this.comparator(t, r) < 0) this.data[e] = r;
			else break;
			e = n;
		}
		return this.data[e] = t, e;
	}
	sink(e) {
		let t = this.data[e], n = this.data.length;
		for (; 2 * e + 1 < n;) {
			let r = 2 * e + 1;
			if (r < n - 1 && this.comparator(this.data[r + 1], this.data[r]) < 0 && r++, this.comparator(t, this.data[r]) <= 0) break;
			this.data[e] = this.data[r], e = r;
		}
		return this.data[e] = t, e;
	}
}, Gc = class extends e {
	constructor(...e) {
		super(...e), this.previousPosition = {
			x: 0,
			y: 0
		}, this.selected = !1, this.collapsed = !1, this.color = "";
	}
	isCollapsed() {
		return this.collapsed;
	}
	setCollapsed(e) {
		this.collapsed = e;
	}
	isSelected() {
		return this.selected;
	}
	getColor() {
		return this.color;
	}
	setSelected(e) {
		this.selected = e;
		for (let t of this.children) t.setSelected(e);
	}
	collapseAll() {
		for (let e of this.children) e.setCollapsed(!0), e.collapseAll();
	}
	collapse() {
		for (let e of this.children) e.setCollapsed(!0);
	}
	expandAll() {
		this.expand(100);
	}
	expand(e) {
		if (e > 0 && this.children.length > 0) for (let t of this.children) t.setCollapsed(!1), t.expand(e - 1);
	}
	setColor(e) {
		this.color = e;
		for (let t of this.children) t.setColor(e);
	}
}, Kc = class e {
	static {
		this.idCounter = 0;
	}
	preprocessData(t) {
		let n = [];
		if (t.children) for (let e of t.children) n.push(this.preprocessData(e));
		return new Gc(t.id || ++e.idCounter, t.name || "", n, t.count, t.selfCount, t.extra);
	}
}, qc = class {
	constructor(e, t, n = new Uc()) {
		this.element = e, this.nodeId = 0, this.zoomScale = 1, this.navigating = !1, this.settings = this.fillOptions(n), this.settings.enableTooltips && (this.tooltip = bc.create(this.element, this.settings.tooltipContainer));
		let r = bo(new Kc().preprocessData(t));
		r.sum((e) => e.children.length > 0 ? 0 : e.count), this.widthScale = Y().range([this.settings.minNodeSize, this.settings.maxNodeSize]), this.treeLayout = Vo().nodeSize([2, 10]).separation((e, t) => {
			if (e.data.isCollapsed() || t.data.isCollapsed()) return 0;
			let n = (this.computeNodeSize(e) + this.computeNodeSize(t)) / 2 + 4;
			return e.parent === t.parent ? n : n + 4;
		}), this.data = this.treeLayout(r).descendants(), this.root = this.data[0], xc.clear(this.element, this.settings.tooltipContainer), this.svg = R(this.element).append("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg").attr("viewBox", `0 0 ${this.settings.width} ${this.settings.height}`).attr("width", this.settings.width).attr("height", this.settings.height).style("font-family", "'Helvetica Neue', Helvetica, Arial, sans-serif"), this.zoomListener = ac().extent([[0, 0], [this.settings.width, this.settings.height]]).scaleExtent([.1, 3]).on("start", () => {
			this.navigating = !0, this.tooltipOut();
		}).on("zoom", (e) => {
			this.zoomScale = e.transform.k, this.visElement.attr("transform", e.transform.toString());
		}).on("end", () => {
			this.navigating = !1;
		}), this.visElement = this.svg.call(this.zoomListener).append("g"), this.render(this.root);
	}
	reset() {
		this.render(this.data[0]);
	}
	fillOptions(e = void 0) {
		let t = new Uc();
		return Object.assign(t, e);
	}
	render(e) {
		this.widthScale.domain([0, e.data.count]), this.root = e, this.root.x = this.settings.height / 2, this.root.y = 0, this.root.data.setSelected(!0);
		let t = (e, n) => {
			if (e.data.setColor(this.settings.colorProvider(e.data, n)), n < this.settings.colorProviderLevels && e.children) for (let r of e.children) t(r, n + 1);
		};
		this.root.children?.forEach((e, n) => {
			t(e, 1);
		}), this.settings.enableExpandOnClick ? (this.root.data.collapseAll(), this.initialExpand(this.root)) : this.root.data.expandAll(), this.update(e), this.centerRoot(e);
	}
	centerRoot(e) {
		let [t, n] = [-e.y, -e.x];
		t = t * this.zoomScale + this.settings.width / 4, n = n * this.zoomScale + this.settings.height / 2, this.visElement.transition().duration(this.settings.animationDuration).attr("transform", `translate(${t},${n})scale(${this.zoomScale})`).on("end", () => this.zoomListener.transform(this.svg, Ys.translate(t, n).scale(this.zoomScale)));
	}
	initialExpand(e) {
		if (!this.settings.enableAutoExpand) {
			e.data.expand(this.settings.levelsToExpand);
			return;
		}
		e.data.expand(1);
		let t = e.data.count * (this.settings.enableAutoExpand ? this.settings.autoExpandValue : .8), n = new Wc([...e.children || []], (e, t) => t.data.count - e.data.count);
		for (; t > 0 && n.size() > 0;) {
			let e = n.remove();
			t -= e.data.count, e.data.expand(1), e.children?.forEach((e, t) => {
				n.add(e);
			});
		}
	}
	update(e) {
		let t = this.treeLayout(this.root), n = t.descendants().reverse().filter((e) => !e.data.isCollapsed()), r = t.links().filter((e) => !e.target.data.isCollapsed() && !e.source.data.isCollapsed());
		n.forEach((e) => e.y = e.depth * this.settings.nodeDistance);
		let i = this.visElement.selectAll("g.node").data(n, (e) => e.data.id || (e.data.id = ++this.nodeId)), a = i.enter().append("g").attr("class", "node").style("cursor", "pointer").attr("transform", `translate(${e.y || 0},${e.data.previousPosition.x || 0})`).on("click", (e, t) => this.click(e, t)).on("mouseover", (e, t) => this.tooltipIn(e, t)).on("mouseout", () => this.tooltipOut()).on("contextmenu", (e, t) => this.rightClick(e, t));
		a.append("circle").attr("r", 1e-6).style("stroke-width", "1.5px").style("stroke", (e) => this.settings.nodeStrokeColor(e.data)).style("fill", (e) => this.settings.nodeFillColor(e.data));
		let o = Y().range([0, 2 * Math.PI]), s = js().innerRadius(0).outerRadius((e) => this.computeNodeSize(e)).startAngle(0).endAngle((e) => o(e.data.selfCount / e.data.count) || 0);
		this.settings.enableInnerArcs && a.append("path").attr("class", "innerArc").attr("d", s).style("fill", (e) => this.settings.nodeStrokeColor(e.data)).style("fill-opacity", 0), this.settings.enableLabels && a.append("text").attr("x", (e) => e.children ? -10 : 10).attr("dy", ".35em").attr("text-anchor", (e) => e.children ? "end" : "start").text((e) => this.settings.getLabel(e.data)).style("font", "10px sans-serif").style("fill-opacity", 1e-6);
		let c = a.merge(i).transition().duration(this.settings.animationDuration).attr("transform", (e) => `translate(${e.y}, ${e.x})`);
		c.select("circle").attr("r", (e) => this.computeNodeSize(e)).style("fill-opacity", (e) => e.children && e.children[0].data.isCollapsed() ? 1 : 0).style("stroke", (e) => this.settings.nodeStrokeColor(e.data)).style("fill", (e) => this.settings.nodeFillColor(e.data)), this.settings.enableInnerArcs && c.select(".innerArc").style("fill-opacity", 1), this.settings.enableLabels && c.select("text").style("fill-opacity", 1);
		let l = i.exit().transition().duration(this.settings.animationDuration).attr("transform", (t) => `translate(${e.y},${e.x})`).remove();
		l.select("circle").attr("r", 1e-6), l.select("path").style("fill-opacity", 1e-6), l.select("text").style("fill-opacity", 1e-6);
		let u = this.visElement.selectAll("path.link").data(r, (e) => e.target.data.id), d = Vs().x((e) => e.y).y((e) => e.x), f = (e, t) => {
			let n = {
				x: e,
				y: t
			};
			return d({
				source: n,
				target: n
			});
		};
		u.enter().insert("path", "g").attr("class", "link").style("fill", "none").style("stroke-opacity", "0.5").style("stroke-linecap", "round").style("stroke", (e) => this.settings.linkStrokeColor(e)).style("stroke-width", 1e-6).attr("d", () => f(e.data.previousPosition.x, e.data.previousPosition.y)).merge(u).transition().duration(this.settings.animationDuration).attr("d", d).style("stroke", this.settings.linkStrokeColor).style("stroke-width", (e) => e.source.data.isSelected() ? this.widthScale(e.target.data.count) + "px" : "4px"), u.exit().transition().duration(this.settings.animationDuration).style("stroke-width", 1e-6).attr("d", () => f(e.x, e.y)).remove(), n.forEach((e) => {
			e.data.previousPosition = {
				x: e.x,
				y: e.y
			};
		});
	}
	computeNodeSize(e) {
		return e.data.isSelected() ? this.widthScale(e.data.count) / 2 : 2;
	}
	click(e, t) {
		this.settings.enableExpandOnClick && (e.defaultPrevented || (e.shiftKey ? t.data.expandAll() : t.children && t.children.some((e) => !e.data.isCollapsed()) ? t.data.collapseAll() : t.data.expand(this.settings.levelsToExpand), this.update(t), this.centerRoot(t)));
	}
	tooltipIn(e, t) {
		if (this.settings.enableTooltips && this.tooltip && !this.navigating) {
			let n = this.settings.getTooltip(t.data);
			this.tooltipTimer = window.setTimeout(() => this.tooltip.show(e, n), 1e3);
		}
	}
	tooltipOut() {
		this.settings.enableTooltips && this.tooltip && (clearTimeout(this.tooltipTimer), this.tooltip.hide());
	}
	rightClick(e, t) {
		this.settings.enableRightClick && this.render(t);
	}
}, Jc = class {
	constructor() {
		this.padding = {
			top: 10,
			right: 10,
			bottom: 10,
			left: 10
		};
	}
}, Yc = class {
	constructor() {
		this.padding = {
			top: 10,
			right: 10,
			bottom: 10,
			left: 10
		}, this.titleFontSize = 24, this.labelFontSize = 16, this.symbolSize = 16, this.columns = 3, this.width = 300, this.rowSpacing = 5, this.columnSpacing = 20;
	}
}, Xc = class extends oc {
	constructor(...e) {
		super(...e), this.height = 800, this.orientation = "vertical", this.barHeight = 75, this.className = "barplot", this.maxItems = 20, this.font = "\"Roboto\", sans-serif", this.displayMode = "relative", this.showBarLabel = !0, this.showValuesInBars = !0, this.barLabelWidth = 150, this.valuesInBarsFontSize = 12, this.chart = new Jc(), this.legend = new Yc(), this.enableTooltips = !0, this.highlightOnHover = !0, this.mouseIn = () => {}, this.mouseMove = () => {}, this.mouseOut = () => {}, this.getTooltip = (e, t, n) => `
            <style>
                .unipept-tooltip {
                    padding: 10px;
                    border-radius: 5px; 
                    background: rgba(0, 0, 0, 0.8); 
                    color: #fff;
                }
                
                .unipept-tooltip div {
                    font-family: "Roboto", sans-serif;
                }
            </style>
            <div class="unipept-tooltip">
                <div style="font-size: 20px; margin-bottom: 8px;">
                    ${this.getTooltipTitle(e, t, n)}
                </div>
                <div>
                    ${this.getTooltipText(e, t, n)}
                </div>
            </div>
        `, this.getTooltipTitle = (e, t, n) => e[t].items[n].label, this.getTooltipText = (e, t, n) => {
			let r = [], i = e[t].items[n].label;
			for (let t of e) {
				let e = t.items.find((e) => e.label === i), n;
				n = e ? this.displayMode === "absolute" ? `${e.counts} hits` : `${e.counts.toFixed(2)}%` : "Not present", r.push(`<span style="font-weight: 600;">${t.label}: </span><span>${n}</span>`);
			}
			return r.map((e) => `<div style="margin-bottom: 4px;">${e}</div>`).join("\n");
		};
	}
}, Zc = class {
	computeMaxItemsInBars(e, t) {
		let n = [...e[0].items].sort((e, t) => t.counts - e.counts);
		t !== void 0 && (n = n.splice(0, t));
		let r = e.map((e) => {
			let t = 0, r = [];
			for (let i of e.items) n.findIndex((e) => e.label === i.label) >= 0 ? r.push(i) : t += i.counts;
			let i = r.sort((e, t) => n.findIndex((t) => t.label === e.label) - n.findIndex((e) => e.label === t.label));
			return {
				label: e.label,
				items: i,
				otherCount: t
			};
		});
		return r.every((e) => e.otherCount === 0) ? r.map(({ label: e, items: t }) => ({
			label: e,
			items: t
		})) : r.map(({ label: e, items: t, otherCount: n }) => ({
			label: e,
			items: [...t, {
				label: "Other",
				counts: n
			}]
		}));
	}
	convertAbsoluteToRelative(e) {
		return e.map((e) => {
			let t = e.items.reduce((e, t) => e + t.counts, 0);
			return {
				label: e.label,
				items: e.items.map((e) => ({
					label: e.label,
					counts: t === 0 ? 0 : e.counts / t * 100
				}))
			};
		});
	}
}, Qc = 40, $c = 5, el = 10, tl = class {
	constructor(e, t, n = new Xc()) {
		this.element = e, this.settings = this.fillOptions(n);
		let r = new Zc();
		this.data = r.computeMaxItemsInBars(t, this.settings.maxItems), this.settings.displayMode === "relative" && (this.data = r.convertAbsoluteToRelative(this.data)), this.settings.enableTooltips && (this.tooltip = bc.create(this.element, this.settings.tooltipContainer)), this.renderBarplot();
	}
	resize(e) {
		this.settings.width = e, this.renderBarplot();
	}
	fillOptions(e = void 0) {
		let t = new Xc();
		return Object.assign(t, e), t.chart = this.fillGroup(new Jc(), e?.chart), t.legend = this.fillGroup(new Yc(), e?.legend), t;
	}
	fillGroup(e, t) {
		let n = Object.assign({}, e.padding, t?.padding);
		return Object.assign(e, t, { padding: n });
	}
	get contentHeight() {
		let e = this.settings.chart.padding.top + this.plotAreaHeight + $c + Qc + this.settings.chart.padding.bottom, t = this.settings.orientation === "horizontal" ? this.legendHeight : this.plotAreaHeight + Qc + this.legendHeight;
		return Math.max(e, t);
	}
	get plotAreaHeight() {
		return this.settings.barHeight * this.data.length;
	}
	get legendHeight() {
		let e = this.settings.legend, t = new Set(this.data.flatMap((e) => e.items.map((e) => e.label))).size, n = Math.ceil(t / e.columns), r = Math.max(e.symbolSize, e.labelFontSize);
		return e.padding.top + e.titleFontSize + el + n * r + Math.max(n - 1, 0) * e.rowSpacing + e.padding.bottom;
	}
	renderBarplot() {
		xc.clear(this.element, this.settings.tooltipContainer);
		let e = R(this.element).append("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg").attr("viewBox", `0 0 ${this.settings.width} ${this.contentHeight}`).attr("width", this.settings.width).attr("height", this.contentHeight).attr("overflow", "hidden").style("font-family", this.settings.font);
		this.initCss();
		let t = this.settings.font, n = this.settings.chart.padding, r = this.settings.orientation == "horizontal", i = this.settings.legend.padding, a = this.settings.legend.width, o = this.settings.legend.titleFontSize, s = this.settings.legend.labelFontSize, c = this.settings.legend.symbolSize, l = this.settings.legend.rowSpacing, u = this.settings.legend.columnSpacing, d = this.settings.legend.columns, f, p, m, h, g, _, v;
		r ? (f = this.settings.width - n.left - n.right - a, m = i.top, p = n.left + f + n.right + i.left, h = Math.max(c, s), g = a - i.left - i.right - c - 10, v = a - i.left - i.right) : (f = this.settings.width - n.left - n.right, m = this.plotAreaHeight + i.top + Qc, p = i.left, h = Math.max(c, s), _ = this.settings.width - i.left - i.right, v = Math.floor((_ - Math.max(d - 1, 0) * u) / d), g = v - c - 10);
		let y = this.settings.barLabelWidth, x = f;
		this.settings.showBarLabel ? x = f - y - 10 : y = 0;
		let S = e.append("g"), C = Ks().keys(Array.from(new Set(this.data.flatMap((e) => e.items.map((e) => e.label))))).value((e, t) => e.items.find((e) => e.label === t)?.counts ?? 0)(this.data), w = Y().domain([0, b(C, (e) => b(e, (e) => e[1])) || 0]).range([0, x]), T = Xo().domain(this.data.map((e, t) => t.toString())).range([0, this.plotAreaHeight]).paddingInner(.1).paddingOuter(0), E = [
			"#9e0142",
			"#c72e4c",
			"#d53e4f",
			"#eb5c48",
			"#f46d43",
			"#fba35b",
			"#fdae61",
			"#fee08b",
			"#ffffbf",
			"#e6f598",
			"#b5e3a5",
			"#8dd380",
			"#66c2a5",
			"#4dacb1",
			"#3288bd",
			"#1f78b4",
			"#5e4fa2",
			"#6a3d9a",
			"#984ea3",
			"#df7ab4"
		];
		this.settings.maxItems && (E[this.settings.maxItems % (this.data[0].items.length + 1)] = "#acaaaa");
		let D = Yo().domain(Array.from(new Set(this.data.flatMap((e) => e.items.map((e) => e.label))))).range(E);
		this.settings.showBarLabel && S.append("g").attr("class", "barLabels").selectAll("text").data(this.data).join("text").attr("x", n.left).attr("y", (e, t) => n.top + (T(t.toString()) || 0) + T.bandwidth() / 2).attr("dy", ".35em").attr("font-family", t).attr("font-size", 18).text((e) => {
			if (e.label.length * (18 * .6) > y) {
				let t = Math.floor(y / (18 * .6));
				return e.label.substring(0, t - 3) + "...";
			}
			return e.label;
		});
		let O = Array.from({ length: this.data.length }, () => []);
		for (let e of C) {
			let t = e.key;
			for (let n = 0; n < e.length; n++) O[n].push({
				barIndex: n,
				title: t,
				shape: e[n]
			});
		}
		S.append("g").selectAll("g").data(O).join("g").selectAll("g").data((e) => e).join((e) => {
			let r = e.append("g");
			return r.append("rect").attr("fill", (e) => D(e.title)).attr("x", (e) => n.left + y + 10 + Math.floor(w(e.shape[0]))).attr("y", (e) => n.top + (T(e.barIndex.toString()) || 0)).attr("width", (e) => Math.floor(w(e.shape[1])) - Math.floor(w(e.shape[0]))).attr("height", Math.floor(T.bandwidth())), this.settings.showValuesInBars && r.append("text").attr("data-key", (e) => e.title).attr("x", (e) => {
				let t = Math.floor(w(e.shape[0])), r = Math.floor(w(e.shape[1]));
				return n.left + y + 10 + t + (r - t) / 2;
			}).attr("y", (e) => n.top + (T(e.barIndex.toString()) || 0) + T.bandwidth() / 2).attr("dy", ".35em").attr("text-anchor", "middle").attr("fill", (e) => {
				let t = jn(D(e.title));
				return (.299 * t.r + .587 * t.g + .114 * t.b) / 255 < .5 ? "white" : "#171717";
			}).attr("font-family", t).attr("font-size", this.settings.valuesInBarsFontSize).attr("font-weight", 600).text((e) => {
				let t = e.shape[1] - e.shape[0];
				return Math.floor(w(e.shape[1])) - Math.floor(w(e.shape[0])) < 30 ? "" : this.settings.displayMode === "relative" ? `${t.toFixed(1)}%` : t;
			}), r;
		}).classed("barplot-item", !0).attr("data-bar-item", (e) => e.title).on("mouseover", (e, t) => {
			let n = this.data[t.barIndex].items.findIndex((e) => e.label === t.title);
			this.mouseIn(e, t.barIndex, n);
		}).on("mousemove", (e, t) => {
			let n = this.data[t.barIndex].items.findIndex((e) => e.label === t.title);
			this.mouseMove(e, t.barIndex, n);
		}).on("mouseout", (e, t) => {
			let n = this.data[t.barIndex].items.findIndex((e) => e.label === t.title);
			this.mouseOut(e, t.barIndex, n);
		}), S.append("g").attr("transform", `translate(${n.left + y + 10}, ${n.top + this.plotAreaHeight + $c})`).call(M(w)).attr("font-size", "12px").append("text").attr("font-family", t).attr("fill", "black").attr("x", x / 2).attr("y", Qc).attr("text-anchor", "middle").attr("font-size", 14).text(this.settings.displayMode === "relative" ? "Percentage" : "Count");
		let k = S.append("g").attr("font-family", t).attr("font-size", s).selectAll("g").data(D.domain()).join("g").classed("legend-item", !0).attr("data-legend-entry", (e) => e).attr("transform", (e, t) => `translate(${t % d * v + Math.max(t % d - 1, 0) * u}, ${Math.floor(t / d) * (h + l) + o + el + m})`);
		S.append("text").attr("font-family", t).attr("font-size", o).attr("dominant-baseline", "hanging").attr("x", p).attr("y", m).text("Legend"), k.append("rect").attr("x", p).attr("width", c).attr("height", c).attr("rx", 5).attr("fill", D), k.append("text").attr("x", p + c + 10).attr("y", s / 2).attr("dy", "0.35em").text((e) => {
			if (e.length * (s * .6) > g) {
				let t = Math.floor(g / (s * .6));
				return e.substring(0, t - 3) + "...";
			}
			return e;
		});
	}
	initCss() {
		let e = this.settings.className;
		Ic.applyStyle(this.element, e, `
.${e} .barplot-item-highlighted {
    opacity: 0.5;
    transition: opacity 0.2s ease-in-out;
    font-size: 20px;
}

.${e} .legend-item-highlighted {
    opacity: 0.5;
    transition: opacity 0.2s ease-in-out;
}
`);
	}
	mouseIn(e, t, n) {
		let r = this.data[t].items[n];
		if (this.settings.mouseIn(this.data, t, n, {
			x: e.clientX,
			y: e.clientY
		}), this.settings.enableTooltips && this.tooltip && this.tooltip.show(e, this.settings.getTooltip(this.data, t, n)), this.settings.highlightOnHover) {
			let e = R(this.element);
			e.selectAll(".barplot-item").classed("barplot-item-highlighted", (e) => e.title !== r.label), e.selectAll(".legend-item").classed("legend-item-highlighted", (e) => e !== r.label);
		}
	}
	mouseMove(e, t, n) {
		this.settings.mouseMove(this.data, t, n, {
			x: e.clientX,
			y: e.clientY
		}), this.settings.enableTooltips && this.tooltip && this.tooltip.move(e);
	}
	mouseOut(e, t, n) {
		if (this.settings.mouseOut(this.data, t, n), this.settings.enableTooltips && this.tooltip && this.tooltip.hide(), this.settings.highlightOnHover) {
			let e = R(this.element);
			e.selectAll(".barplot-item").classed("barplot-item-highlighted", !1), e.selectAll(".legend-item").classed("legend-item-highlighted", !1);
		}
	}
};
//#endregion
export { tl as Barplot, Jc as BarplotChartSettings, Yc as BarplotLegendSettings, Xc as BarplotSettings, jc as ColorPalette, Lc as ColorUtils, e as DataNode, dc as EuclidianDistanceMetric, kc as Heatmap, pc as HeatmapLegendSettings, mc as HeatmapSettings, fc as MoloReorderer, Ac as PearsonCorrelationMetric, Mc as StringUtils, Rc as Sunburst, Nc as SunburstSettings, sc as Transition, cc as TreeNode, Vc as Treemap, zc as TreemapSettings, qc as Treeview, Uc as TreeviewSettings, uc as UPGMAClusterer };

//# sourceMappingURL=unipept-visualizations.js.map