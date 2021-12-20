window.addEventListener('load', main, false);
function main() {

    var art = [[0, -82.6,  25.7,  157.3,  14.9,  9.80, -127.9,    9.5,   20.1,     251.6,    1541.4,     146.2],
    [24, -20.5, 146.7,   12.9, 110.7, 33.68,  -19.2,    4.7,   11.5,    4942.7,     433.2,    3729.2],
    [48, -18.9, 144.3,   12.2, 123.2, 49.37,  -18.1,    4.9,   12.8,    7126.8,     602.6,    6082.6],
    [72, -18.3, 143.4,   12.0, 128.0, 60.15,  -17.7,    5.0,   13.2,    8625.9,     719.8,    7698.8],
    [96, -18.0, 142.9,   11.8, 130.4 ,67.65,  -17.5,    5.0,   13.4,    9669.0,     801.5,    8823.7],
   [120, -17.9, 142.7,   11.8, 131.8, 72.90,  -17.4,    5.0,   13.5,   10399.0,     858.8,    9611.0],
   [144, -17.8, 142.5,   11.7, 132.7, 76.58,  -17.3,    5.0,   13.6,   10911.3,     899.0,   10163.4],
   [168, -17.7, 142.4,   11.7, 133.3, 79.17,  -17.2,    5.0,   13.7,   11271.3,     927.2,   10551.7],
   [192, -17.6, 142.3,   11.7, 133.7, 80.99,  -17.2,    5.0,   13.7,   11524.6,     947.1,   10824.9],
   [216, -17.6, 142.2,   11.7, 133.9, 82.27,  -17.2,    5.1,   13.7,   11702.9,     961.1,   11017.2],
   [240, -17.6, 142.2,   11.7, 134.1, 83.17,  -17.2,    5.1,   13.7,   11828.4,     971.0,   11152.6]];
    var t_true = [0, 24, 48, 72, 96, 120, 144, 168, 192, 216, 240];
    var na_true = [25.7, 146.7, 144.3, 143.4, 142.9, 142.7, 142.5, 142.4, 142.3, 142.2, 142.2,];
    var k_true = [157.3, 12.9, 12.2, 12.0, 11.8, 11.8, 11.7, 11.7, 11.7, 11.7, 11.7,];
    var cl_true = [14.9, 110.7, 123.2, 128.0, 130.4, 131.8, 132.7, 133.3, 133.7, 133.9, 134.1,];

    var t_max = 300;
    var hp = 1;
    var lt;
    var indic;
    var mt = [];
    var mu = [];
    var mk = [];
    var mna = [];
    var mcl = [];
    var ml = [];
    var mva = [];
    var mz = [];
    var na2;
    var k2;
    var cl2;
    var V2;
    var beta2;
    var l2;
    var kp;
    var tind;
    var naind;
    var mast;
    var B0;
    var osor;
    var pp;
    var hi;

    // READING DATA
    var na0 = 140.0;
    var k0 = 5.0;
    var cl0 = 155.0;
    B0 = 0.0;
    kv = 1.000;
    na1 = 25.7;
    k1 = 157.3;
    cl1 = 14.9;
    beta = 0.050;
    gamma = 1.50;
    pna = 0.05;
    pk = 3.33;
    pcl = 0.035;
    inc = 0.001;
    ikc = 0;
    inkcc = 0;
    hp = 1;
    kb = 0;

    // idk assignment
    var na = na1;
    var k = k1;
    var cl = cl1;
    var BB = B0;
    var A1 = 1;
    na = na*kv;
    k = k*kv;
    cl = cl*kv;
    mast = (k0+na0+cl0+BB-k-na-cl);
    var V = A1/mast;
    z = V*(cl-na-k)/A1;
    na2 = na;
    k2 = k;
    cl2 = cl;
    beta2 = beta;
    V2 = V;
    var prna;
    var prk;
    var prcl;
    var h = 0.1;
    var t = 0;
    lt =  1;
    
    var unach = -5.5;
    var ukon = 0.2;

    if ((na2<=0) || (k2<=0) || (cl2<=0) || (V2<=0) || (V<=0) || (z>=0)) {
        throw 'BAD INITIAL DATA';
    }

    function F(u, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma) {
        return Math.exp(u) - (pna*na0 + pk*k0 + pcl*cl + beta*na*(1-1/gamma)/u)/(pna*na + pk*k + pcl*cl0 + beta*na*(1-1/gamma)/u);
    }

    function solveU(na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma) {
        var eps = 0.0001;
        var a = unach;
        var b = -eps;
        while (F(a, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma)*F(b, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma) > 0) {
            b-=eps;
            if (a>b) {
                throw 'Can\'t solve F(u)=0';
            }
        }
        while ((b-a)/2 > eps) {
            c = (b+a)/2;
            if (F(c, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma)*F(b, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma)<0) {
                a = c;
            } else {
                b = c;
            }
        }
        return c;
    }

    function calc() {
        var u = solveU(na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma);
        mt.push(t);
        mu.push(u);
        mna.push(na);
        mk.push(k);
        mcl.push(cl);
        mva.push(V);
        mz.push(z);

        var i=0;
        do {
            // solve eq for u    vichu in pascal
            var u = solveU(na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma);
            //console.log(u);
            // find devariatives for Na, K, Cl   PROIZBEZ in pascal
            var g = 1-Math.exp(u);
            if ((cl0!=0) && (k0!=0) && (na0!=0)) {
                var jnc = inc*(na0*cl0-na*cl);
                var jkc = ikc*(k0*cl0-k*cl);
                var jnkcc = inkcc*(na0*k0*cl0*cl0-na*k*cl*cl); 
            }
            V = A1/mast;
            var fna = pna*u*(na*Math.exp(u)-na0)/g + jnc + jnkcc - beta*na;
            var fk = pk*u*(k*Math.exp(u)-k0)/g + jnkcc + jkc + beta/gamma*na;
            var fcl = pcl*u*(cl-cl0*Math.exp(u))/g + jnc + 2*jnkcc + jkc;
            mast = (k0+na0+cl0+BB-k-na-cl);
            z = V*(cl-na-k)/A1;
            prna = fna*(1-na/(na0+k0+cl0+BB)) - na/(na0+k0+cl0+BB)*(fk+fcl);
            prk = fk*(1-k/(na0+k0+cl0+BB)) - k/(na0+k0+cl0+BB)*(fna+fcl);
            prcl = fcl*(1-cl/(na0+k0+cl0+BB)) - cl/(na0+k0+cl0+BB)*(fk+fna);
            na = na+h*prna;
            if ((na<=0.1) && (indic==0)) {
                indic = 1;
                tind = t;
                naind = na;
                throw `LOW SODIUM t=${tind.toFixed(4)} na=${naind.toFixed(4)}`;
            }
            k = k+prk*h;
            cl = cl+prcl*h;
            mast = (k0+na0+cl0+BB-k-na-cl);
            z = (cl-na-k)/mast;
            V = A1/mast;
            t = t+h;
            beta = beta-kb*h;

            

            if (i%hp==0) {
                mt.push(t);
                mu.push(u);
                mna.push(na);
                mk.push(k);
                mcl.push(cl);
                mva.push(V);
                mz.push(z);
                mast = (k0+na0+cl0+BB-k-na-cl);
                //console.log(fna+fk-fcl);
                //console.log(t.toFixed(2), (u*26.7).toFixed(2), na.toFixed(2), k.toFixed(2), cl.toFixed(2), (V*1000).toFixed(6), prna.toFixed(3), prk.toFixed(3), prcl.toFixed(3));
                lt++;
            }
            i++;
        } while ((Math.abs(prna)+Math.abs(prk)+Math.abs(prcl)>1e-8) && (t<t_max))
        Plotly.newPlot('tester', [
            {
                x: mt,
                y: mna,
                mode: 'lines',
                name: 'Na',
                line: {
                    color: 'red',
                    width: 2,
                }
            },
            {
                x: mt,
                y: mk,
                mode: 'lines',
                name: 'K',
                line: {
                    color: 'blue',
                    width: 2,
                }
            },
            {
                x: mt,
                y: mcl,
                mode: 'lines',
                name: 'Cl',
                line: {
                    color: 'green',
                    width: 2,
                }
            },
            /*{
                x: t_true,
                y: na_true,
                mode: 'markers',
                name: 'Na pascal',
                marker: {
                    color: 'red',
                    size: 6,
                }
            },
            {
                x: t_true,
                y: k_true,
                mode: 'markers',
                name: 'K pascal',
                marker: {
                    color: 'blue',
                    size: 6,
                }
            },
            {
                x: t_true,
                y: cl_true,
                mode: 'markers',
                name: 'Cl pascal',
                marker: {
                    color: 'green',
                    size: 6,
                }
            },*/
        ]);
    }

    calc()
        


    function plot(X, Na, K, Cl, width=900, height=400) {
        let canvas = document.createElement('canvas');
        canvas.height = height;
        canvas.width = width;
        canvas.style = 'border: 1px solid black;';
        document.body.appendChild(canvas);
        let ctx = canvas.getContext('2d');
        let X_max = Math.max(...X);
        let X_min = Math.min(...X);
        let Y_max = 200//Math.max(...Na, ...K, ...Cl)*1.1;
        let Y_min = 0//Math.min(...Na, ...K, ...Cl)-Y_max*0.1;
        let X_scale = width/(X_max-X_min);
        let Y_scale = height/(Y_max-Y_min);
        ctx.beginPath();
        ctx.moveTo(0, height+Y_min*Y_scale);
        ctx.lineTo(width, height+Y_min*Y_scale);
        ctx.stroke();
        console.log(`X limits (${X_min}; ${X_max})\nY limits (${Y_min}; ${Y_max})`);
        ctx.lineWidth = 3;
        ctx.beginPath();
        ctx.moveTo(X[0], height-(Na[0]-Y_min)*Y_scale);
        for (let i=0; i<X.length; i++) {
            if (X[i]!=NaN || Na[i]!=NaN) {
                ctx.lineTo((X[i]-X_min)*X_scale, height-(Na[i]-Y_min)*Y_scale);
            }
        }
        ctx.strokeStyle = "red";
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(X[0], height-(K[0]-Y_min)*Y_scale);
        for (let i=0; i<X.length; i++) {
            if (X[i]!=NaN || K[i]!=NaN) {
                ctx.lineTo((X[i]-X_min)*X_scale, height-(K[i]-Y_min)*Y_scale);
            }
        }
        ctx.strokeStyle = "blue";
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(X[0], height-(Cl[0]-Y_min)*Y_scale);
        for (let i=0; i<X.length; i++) {
            if (X[i]!=NaN || Cl[i]!=NaN) {
                ctx.lineTo((X[i]-X_min)*X_scale, height-(Cl[i]-Y_min)*Y_scale);
            }
        }
        ctx.strokeStyle = "green";
        ctx.stroke();
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(art[0][0], height-(Na[0]-Y_min)*Y_scale);
        for (let i=0; i<art.length; i++) {
            ctx.lineTo((art[i][0]-X_min)*X_scale, height-(art[i][2]-Y_min)*Y_scale);
            ctx.arc((art[i][0]-X_min)*X_scale, height-(art[i][2]-Y_min)*Y_scale, 3, 0, 2*Math.PI);
        }
        ctx.strokeStyle = "red";
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(art[0][0], height-(K[0]-Y_min)*Y_scale);
        for (let i=0; i<art.length; i++) {
            ctx.lineTo((art[i][0]-X_min)*X_scale, height-(art[i][3]-Y_min)*Y_scale);
            ctx.arc((art[i][0]-X_min)*X_scale, height-(art[i][3]-Y_min)*Y_scale, 3, 0, 2*Math.PI);
        }
        ctx.strokeStyle = "blue";
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(art[0][0], height-(Cl[0]-Y_min)*Y_scale);
        for (let i=0; i<art.length; i++) {
            ctx.lineTo((art[i][0]-X_min)*X_scale, height-(art[i][4]-Y_min)*Y_scale);
            ctx.arc((art[i][0]-X_min)*X_scale, height-(art[i][4]-Y_min)*Y_scale, 3, 0, 2*Math.PI);
        }
        ctx.strokeStyle = "green";
        ctx.stroke();
    }

    //plot(mt, mna, mk, mcl);
    

    form_submit.onclick = function (e) {
        e.preventDefault();
        t_max = 300;
        hp = 1;
        mt = [];
        mu = [];
        mk = [];
        mna = [];
        mcl = [];
        ml = [];
        mva = [];
        mz = [];

        // READING DATA
        na0 = parseFloat(input_na0.value);
        k0 = parseFloat(input_k0.value);
        cl0 = parseFloat(input_cl0.value);
        B0 = parseFloat(input_b0.value);
        kv = parseFloat(input_kv.value);
        na1 = parseFloat(input_na.value);
        k1 = parseFloat(input_k.value);
        cl1 = parseFloat(input_cl.value);
        beta = parseFloat(input_beta.value);
        gamma = parseFloat(input_gamma.value);
        pna = parseFloat(input_pna.value);
        pk = parseFloat(input_pk.value);
        pcl = parseFloat(input_pcl.value);
        inc = parseFloat(input_inc.value);
        ikc = parseFloat(input_ikc.value);
        inkcc = parseFloat(input_inkcc.value);
        hp = 1;
        kb = 0;

        // idk assignment
        na = na1;
        k = k1;
        cl = cl1;
        BB = B0;
        A1 = 1;
        na = na*kv;
        k = k*kv;
        cl = cl*kv;
        mast = (k0+na0+cl0+BB-k-na-cl);
        V = A1/mast;
        z = V*(cl-na-k)/A1;
        na2 = na;
        k2 = k;
        cl2 = cl;
        beta2 = beta;
        V2 = V;
        h = 0.1;
        t = 0;
        lt =  1;
        calc();
    }
}