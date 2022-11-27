window.addEventListener('load', main, false);
function main() {

    const isNumber = value => typeof value === 'number' && value === value && value !== Infinity && value !== -Infinity;

    // all variables
    var na, k, cl, na0, k0, cl0, t_max, hp, lt, indic, mt, mu, mk, mna, mcl, ml, mva, mz, 
        naC, kC, clC, mosor, mun, muk, mucl, na2, k2, cl2, V2, beta2, tind, naind, mast, 
        B0, osor, kv, na1, k1, cl1, beta, gamma, pna, pk, pcl, inc, ikc, inkcc, kb, BB, A1,
        t, lt, unach, ukon, g, u, prna, prcl, prk, flux, error;

    // initialization
    function init() {
        // read data from form
        na1 = parseFloat(input_na.value);
        k1 = parseFloat(input_k.value);
        cl1 = parseFloat(input_cl.value);
        na0 = parseFloat(input_na0.value);
        k0 = parseFloat(input_k0.value);
        cl0 = parseFloat(input_cl0.value);
        B0 = parseFloat(input_b0.value);
        beta = parseFloat(input_beta.value);
        gamma = parseFloat(input_gamma.value);
        pna = parseFloat(input_pna.value);
        pk = parseFloat(input_pk.value);
        pcl = parseFloat(input_pcl.value);
        inc = parseFloat(input_inc.value);
        ikc = parseFloat(input_ikc.value);
        inkcc = parseFloat(input_inkcc.value);
        kv = parseFloat(input_kv.value);
        hp = parseFloat(input_hp.value);
        t_max = parseFloat(input_tmax.value); 

        // check readed data
        if (!(isNumber(na0) && isNumber(k0) && isNumber(cl0) && isNumber(B0) && isNumber(kv) && isNumber(na1) && isNumber(k1) && isNumber(cl1) && isNumber(beta) && isNumber(gamma) && isNumber(pna) && isNumber(pk) && isNumber(pcl) && isNumber(inc) && isNumber(ikc) && isNumber(inkcc) && isNumber(t_max))) {
            warning_div.style.display = 'inline-block';
            throw 'BAD DATA'
        }

        // zeroing arrays
        mt = [];
        mu = [];
        mk = [];
        mna = [];
        mcl = [];
        ml = [];
        mva = [];
        mz = [];
        naC = [];
        kC = [];
        clC = [];
        mosor = [];
        mun = [];
        muk = [];
        mucl = [];
        flux = {
            'Net' : {
              'Pump': {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'Channel': {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'NC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'KC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'NKCC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
            }, 
            'In' : {
              'Pump': {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'Channel': {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'NC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'KC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'NKCC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
            },
            'Ef' : {
              'Pump': {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'Channel': {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'NC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'KC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
              'NKCC' : {
                'Na': [],
                'K' : [],
                'Cl': [],
              },
            }
          };

        // assertion
        kb = 0;
        na = na1;
        k = k1;
        BB = B0;
        A1 = 1;
        na = na*kv;
        k = k*kv;
        cl = cl1*kv;
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
        lt = 1;
        unach = -5.5;

        naC = [];
        kC = [];
        clC = [];

        if ((na2<=0) || (k2<=0) || (cl2<=0) || (V2<=0) || (V<=0) || (z>=0)) {
            throw 'BAD INITIAL DATA';
        }
    
    }

    function F(u, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma) {
        return Math.exp(u) - (pna*na0 + pk*k0 + pcl*cl + beta*na*(1-1/gamma)/u)/(pna*na + pk*k + pcl*cl0 + beta*na*(1-1/gamma)/u);
    }

    function solveU(na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma, unach, ukon, eps=1e-6) {
        var a = unach;
        var b = ukon;
        while (F(a, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma)*F(b, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma) > 0) {
            b-=eps;
            if (a>b) {
                throw 'Can\'t solve F(u)=0';
            }
        }
        var c = (b+a)/2;
        while ((b-a)/2 > eps) {
            c = (b+a)/2;
            if (F(c, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma)*F(b, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma)<0) {
                a = c;
            } else {
                b = c;
            }
        }
		console.log(F(c, na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma))
        return c;
    }

    function calc_fluxes() {
        var jnc=inc*(na0*cl0-na*cl);
        var jkc=ikc*(k0*cl0-k*cl);
        var jnkcc=inkcc*(na0*k0*cl0*cl0-na*k*cl*cl);
  
        var jnce=inc*(-na*cl);
        var jkce=ikc*(-k*cl);
        var jnkcce=inkcc*(-na*k*cl*cl);
        flux['Net']['Pump']['Na'].push(-beta*na);
        flux['Net']['Channel']['Na'].push(pna*u*(na*Math.exp(u)-na0)/g);
        flux['Net']['NC']['Na'].push(jnc);
        flux['Net']['KC']['Na'].push(0);
        flux['Net']['NKCC']['Na'].push(jnkcc);
        flux['Net']['Pump']['K'].push(beta/gamma*na);
        flux['Net']['Channel']['K'].push(pk*u*(k*Math.exp(u)-k0)/g);
        flux['Net']['NC']['K'].push(0);
        flux['Net']['KC']['K'].push(jkc);
        flux['Net']['NKCC']['K'].push(jnkcc);
        flux['Net']['Pump']['Cl'].push(0);
        flux['Net']['Channel']['Cl'].push(pcl*u*(cl-cl0*Math.exp(u))/g);
        flux['Net']['NC']['Cl'].push(jnc);
        flux['Net']['KC']['Cl'].push(jkc);
        flux['Net']['NKCC']['Cl'].push(2*jnkcc);
        flux['In']['Pump']['Na'].push(0);
        flux['In']['Channel']['Na'].push(-pna*u*na0/g);
        flux['In']['NC']['Na'].push(inc*na0*cl0);
        flux['In']['KC']['Na'].push(0);
        flux['In']['NKCC']['Na'].push(inkcc*na0*k0*cl0*cl0);
        flux['In']['Pump']['K'].push(beta/gamma*na);
        flux['In']['Channel']['K'].push(pk*u*(-k0)/g);
        flux['In']['NC']['K'].push(0);
        flux['In']['KC']['K'].push(ikc*k0*cl0);
        flux['In']['NKCC']['K'].push(inkcc*na0*k0*cl0*cl0);
        flux['In']['Pump']['Cl'].push(0);
        flux['In']['Channel']['Cl'].push(pcl*u*(-cl0*Math.exp(u))/g);
        flux['In']['NC']['Cl'].push(inc*na0*cl0);
        flux['In']['KC']['Cl'].push(ikc*k0*cl0);
        flux['In']['NKCC']['Cl'].push(2*inkcc*na0*k0*cl0*cl0);
        flux['Ef']['Pump']['Na'].push(-beta*na);
        flux['Ef']['Channel']['Na'].push(pna*u*na*Math.exp(u)/g);
        flux['Ef']['NC']['Na'].push(jnce);
        flux['Ef']['KC']['Na'].push(0);
        flux['Ef']['NKCC']['Na'].push(jnkcce);
        flux['Ef']['Pump']['K'].push(0);
        flux['Ef']['Channel']['K'].push(pk*u*k*Math.exp(u)/g);
        flux['Ef']['NC']['K'].push(0);
        flux['Ef']['KC']['K'].push(jkce);
        flux['Ef']['NKCC']['K'].push(jnkcce);
        flux['Ef']['Pump']['Cl'].push(0);
        flux['Ef']['Channel']['Cl'].push(pcl*u*cl/g);
        flux['Ef']['NC']['Cl'].push(jnce);
        flux['Ef']['KC']['Cl'].push(jkce);
        flux['Ef']['NKCC']['Cl'].push(2*jnkcce);
    }

    function DIAPBEZL() {
        let unach1 = unach;
        let ukon1 = ukon;
        let deltau = 0.013;
        let uu = unach1;
        let ff1 = 0.1;
        let ff2 = 0.1;
        do {
            let ch1 = (pcl*cl+pk*k0+pna*na0+beta*na*(1-1/gamma)/uu);
            let zn1 = (pcl*cl0+pk*k+pna*na+beta*na*(1-1/gamma)/uu);
            ff1 = zn1*Math.exp(uu)-ch1;
            uu += deltau;
            ch1 = (pcl*cl+pk*k0+pna*na0+beta*na*(1-1/gamma)/uu);
            zn1 = (pcl*cl0+pk*k+pna*na+beta*na*(1-1/gamma)/uu);
            ff2 = zn1*Math.exp(uu)-ch1;
            if (uu>=ukon1) {
                console.error('RANGE NOT FOUND');
                return [uu-deltau, uu];
            }
        } while (!((ff1*ff2<0) && (ff1!=0) && (ff2!=0)))
		return [uu-deltau, uu];
    }

    function ZEROIN(ax, bx, tol, F) {
		let eps = 1.0;
		let fierstcall = true;
		let a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s;
		if (fierstcall) {
			do {
				eps = eps/2;
				tol1 = 1+eps;
			} while (tol1 <= 1)
		}
    }

    function calculate() {
        unach = -5.5;
        ukon = 0.2;
        error = 1e-10;
        [unach, ukon] = DIAPBEZL();
		console.log(unach, ukon);
        //u = ZEROIN(unach, ukon, error);
        u = solveU(na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma, unach, ukon);
        mt.push(t);
        mu.push(u*26.7);
        mna.push(na);
        mk.push(k);
        mcl.push(cl);
        mva.push(V);
        mz.push(z);
        naC.push(na*V*1000);
        kC.push(k*V*1000);
        clC.push(cl*V*1000);
        g = 1-Math.exp(u);
        osor = (beta/gamma*na)/(ikc*k0*cl0+inkcc*na0*k0*cl0*cl0-pk*u*k0/g);
        mosor.push(osor);
        calc_fluxes();

        if (na!=0)
            mun.push(26.7*(u-Math.log(na0/na)));
        else
            mun.push(0);
        if (k!=0)
            muk.push(26.7*(u-Math.log(k0/k)));
        else 
            muk.push(0);
        if (cl!=0)
            mucl.push(26.7*(-u-Math.log(cl0/cl)));
        else 
            mucl.push(0);

        var i=1;
        do {
			unach = -5.5;
			ukon = 0.2;
			[unach, ukon] = DIAPBEZL();
			console.log(unach, ukon);
            u = solveU(na, k, cl, na0, k0, cl0, pna, pk, pcl, beta, gamma, unach, ukon);
            g = 1-Math.exp(u);
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
            osor = (beta/gamma*na)/(ikc*k0*cl0+inkcc*na0*k0*cl0*cl0-pk*u*k0/g );
            t = t+h;
            beta = beta-kb*h;
            if (i%hp==0) {
                mt.push(t);
                mu.push(u*26.7);
                mna.push(na);
                mk.push(k);
                mcl.push(cl);
                mva.push(V);
                mz.push(z);
                naC.push(na*V*1000);
                kC.push(k*V*1000);
                clC.push(cl*V*1000);
                mosor.push(osor);
                mun.push(26.7*(u-Math.log(na0/na)));
                muk.push(26.7*(u-Math.log(k0/k)));
                mucl.push(26.7*(-u-Math.log(cl0/cl)));
                mast = (k0+na0+cl0+BB-k-na-cl);
                calc_fluxes();
                lt++;
            }
            i++;
        } while ((Math.abs(prna)+Math.abs(prk)+Math.abs(prcl)>1e-8) && (t<t_max))

        Plotly.newPlot('concentration', [
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
        ],
        {
            title: {
              text:'Ion concentration',
              font: {
                family: 'Courier New, monospace',
                size: 24
              },
              xref: 'paper',
              x: 0.05,
            },
            xaxis: {
              title: {
                text: 'Time (min)',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              },
            },
            yaxis: {
              title: {
                text: 'mM',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              }
            }
        });
        
        Plotly.newPlot('content', [
            {
                x: mt,
                y: naC,
                mode: 'lines',
                name: 'Na',
                line: {
                    color: 'red',
                    width: 2,
                }
            },
            {
                x: mt,
                y: kC,
                mode: 'lines',
                name: 'K',
                line: {
                    color: 'blue',
                    width: 2,
                }
            },
            {
                x: mt,
                y: clC,
                mode: 'lines',
                name: 'Cl',
                line: {
                    color: 'green',
                    width: 2,
                }
            },
        ],
        {
            title: {
              text:'Ion content',
              font: {
                family: 'Courier New, monospace',
                size: 24
              },
              xref: 'paper',
              x: 0.05,
            },
            xaxis: {
              title: {
                text: 'Time (min)',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              },
            },
            yaxis: {
              title: {
                text: 'mmol/mol A',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              }
            }
        });

        Plotly.newPlot('volume', [
            {
                x: mt,
                y: mva,
                mode: 'lines',
                name: 'V/A',
                line: {
                    color: 'black',
                    width: 2,
                }
            },
        ],
        {
            title: {
              text:'V/A',
              font: {
                family: 'Courier New, monospace',
                size: 24
              },
              xref: 'paper',
              x: 0.05,
            },
            xaxis: {
              title: {
                text: 'Time (min)',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              },
            },
            yaxis: {
              title: {
                text: 'ml/mmol A',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              }
            }
        });

        Plotly.newPlot('U', [
            {
                x: mt,
                y: mu,
                mode: 'lines',
                name: 'U',
                line: {
                    color: 'black',
                    width: 2,
                }
            },
            {
                x: mt,
                y: mun,
                mode: 'lines',
                name: 'mun',
                line: {
                    color: 'red',
                    width: 2,
                }
            },
            {
                x: mt,
                y: muk,
                mode: 'lines',
                name: 'muk',
                line: {
                    color: 'blue',
                    width: 2,
                }
            },
            {
                x: mt,
                y: mucl,
                mode: 'lines',
                name: 'mucl',
                line: {
                    color: 'green',
                    width: 2,
                }
            },
        ],
        {
            title: {
              text:'U',
              font: {
                family: 'Courier New, monospace',
                size: 24
              },
              xref: 'paper',
              x: 0.05,
            },
            xaxis: {
              title: {
                text: 'Time (min)',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              },
            },
            yaxis: {
              title: {
                text: 'mV',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              }
            }
        });
        Plotly.newPlot('osor', [
            {
                x: mt,
                y: mosor,
                mode: 'lines',
                name: 'OSOR',
                line: {
                    color: 'black',
                    width: 2,
                }
            },
        ],
        {
            title: {
              text:'OSOR',
              font: {
                family: 'Courier New, monospace',
                size: 24
              },
              xref: 'paper',
              x: 0.05,
            },
            xaxis: {
              title: {
                text: 'Time (min)',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              },
            },
            yaxis: {
              title: {
                text: '',
                font: {
                  family: 'Courier New, monospace',
                  size: 18,
                  color: '#7f7f7f'
                }
              }
            }
        });
    }

    form_submit.onclick = function (e) {
        e.preventDefault();
        init();
        calculate();
    }

    function download(filename, text) {
        var element = document.createElement('a');
        element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
        element.setAttribute('download', filename);
      
        element.style.display = 'none';
        document.body.appendChild(element);
      
        element.click();
      
        document.body.removeChild(element);
    }

    download_result.onclick = function(e) {
        var data = {
            't': mt.map((value)=>value.toFixed(1)),
            'U': mu.map((value)=>value.toFixed(1)),
            'na': mna.map((value)=>value.toFixed(1)),
            'k': mk.map((value)=>value.toFixed(1)),
            'cl': mcl.map((value)=>value.toFixed(1)),
            'V/A': mva.map((value)=>value.toFixed(5)),
            'mun': mun.map((value)=>value.toFixed(1)),
            'muk': muk.map((value)=>value.toFixed(1)),
            'mucl': mucl.map((value)=>value.toFixed(1)),
            'naC': naC.map((value)=>value.toFixed(1)),
            'kC': kC.map((value)=>value.toFixed(1)),
            'clC': clC.map((value)=>value.toFixed(1)),
            'z': mz.map((value)=>value.toFixed(2)),
            'OSOR': mosor.map((value)=>value.toFixed(2)), 
        };
        for (const fluxtype of Object.keys(flux)) {
          for (const transfertype of Object.keys(flux[fluxtype])) {
            for (const atomtype of Object.keys(flux[fluxtype][transfertype])) {
              data[`${fluxtype}flux_${transfertype}_${atomtype}`] = flux[fluxtype][transfertype][atomtype].map((value)=>value.toFixed(4));
            }
          }
        }
        var fields = Object.keys(data);
        var result_text = fields.join(' ')+'\n\n';
        for (var i=0; i<data['t'].length; i++) {
            for (const key of fields) {
                result_text += data[key][i]+' '
            }
            result_text = result_text.substring(0, result_text.length-1)+'\n'
        }
        var jnc=inc*(na0*cl0-na*cl);
        var jkc=ikc*(k0*cl0-k*cl);
        var jnkcc=inkcc*(na0*k0*cl0*cl0-na*k*cl*cl);

        var jnce=inc*(-na*cl);
        var jkce=ikc*(-k*cl);
        var jnkcce=inkcc*(-na*k*cl*cl);
        result_text += `\n* na0 k0 cl0 B0 kv na k cl beta gamma\n* ${na0} ${k0} ${cl0} ${B0} ${kv} ${na1} ${k1} ${cl1} ${beta} ${gamma}\n* pna pk pcl inc ikc inkcc hp kb\n* ${pna} ${pk} ${pcl} ${inc} ${ikc} ${inkcc} ${hp} ${kb}\n\n* Net_flux PUMP Channel NC KC NKCC\n* Na ${(-beta*na).toFixed(4)} ${(pna*u*(na*Math.exp(u)-na0)/g).toFixed(4)} ${(jnc).toFixed(4)} ${0} ${(jnkcc.toFixed(4))}\n* K ${(beta/gamma*na).toFixed(4)} ${(pk*u*(k*Math.exp(u)-k0)/g).toFixed(4)} ${0} ${jkc.toFixed(4)} ${jnkcc.toFixed(4)}\n* Cl ${0} ${(pcl*u*(cl-cl0*Math.exp(u))/g).toFixed(4)} ${jnc.toFixed(4)} ${jkc.toFixed(4)} ${(2*jnkcc).toFixed(4)}\n* Influx PUMP IChannel INC IKC INKCC\n* Na ${0} ${(-pna*u*na0/g).toFixed(4)} ${(inc*na0*cl0).toFixed(4)} ${0} ${(inkcc*na0*k0*cl0*cl0).toFixed(4)}\n* K ${(beta/gamma*na).toFixed(4)} ${(pk*u*(-k0)/g).toFixed(4)} ${0} ${(ikc*k0*cl0).toFixed(4)} ${(inkcc*na0*k0*cl0*cl0).toFixed(4)}\n* Cl ${0} ${(pcl*u*(-cl0*Math.exp(u))/g).toFixed(4)} ${(inc*na0*cl0).toFixed(4)} ${(ikc*k0*cl0).toFixed(4)} ${(2*inkcc*na0*k0*cl0*cl0).toFixed(4)}\n* Efflux PUMP EChannel ENC EKC ENKCC\n* Na ${(-beta*na).toFixed(4)} ${(pna*u*na*Math.exp(u)/g).toFixed(4)} ${jnce.toFixed(4)} ${0} ${jnkcce.toFixed(4)}\n* K ${0} ${(pk*u*k*Math.exp(u)/g).toFixed(4)} ${0} ${jkce.toFixed(4)} ${jnkcce.toFixed(4)}\n* Cl ${0} ${(pcl*u*cl/g).toFixed(4)} ${jnce.toFixed(4)} ${jkce.toFixed(4)} ${(2*jnkcce).toFixed(4)}\n\n* z OSOR (A/V)*1000\n* ${z.toFixed(2)} ${((beta/gamma*na)/(ikc*k0*cl0+inkcc*na0*k0*cl0*cl0-pk*u*k0/g )).toFixed(2)} ${((k0+na0+cl0+BB-k-na-cl)).toFixed(2)}`
        download('result.txt', result_text);
    }



    init();
    calculate();


}
