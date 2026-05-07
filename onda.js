const canvas = document.getElementById('acusticSound');
const ctx = canvas.getContext('2d');


const gridR = 800; // Alcance (Eixo X)
const gridZ = 400; // Profundidade (Eixo Y)
canvas.width = gridR;
canvas.height = gridZ;

/**
 * Implementação da Equação paraxial
 * Resolvendo: iu_x + b*u_zz + V*u = 0
 */

const b = 5;
const Lz = 500;
const Lx = 3000;
const dz = 1;
const dx = 5;

/*
Definindo o tamanho do passo:

z = 0:dz:Lz;
x = 0:dx:Lx;

Nz = length(z);
Nx = length(x);
*/

const nz = Math.floor(Lz / dz) + 1;
const nx = Math.floor(Lx / dx) + 1;
const N = nz - 2;

// Array p/ guardar a solução estática de cada momento |u|
const staticSolution = [];

// Auxiliares de números complexos
const C = {
    add: (a, b) => ({ re: a.re + b.re, im: a.im + b.im }),
    sub: (a, b) => ({ re: a.re - b.re, im: a.im - b.im }),
    mul: (a, b) => ({ re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }),
    scale: (a, s) => ({ re: a.re * s, im: a.im * s }),
    mag: (a) => Math.sqrt(a.re * a.re + a.im * a.im),
    div: (a, b) => {
        const den = b.re * b.re + b.im * b.im;
        return { re: (a.re * b.re + a.im * b.im) / den, im: (a.im * b.re - a.re * b.im) / den };
    }
};

// ==========================================
//  Algoritmo de Thomas (Complexo)
// ==========================================
function solveThomas(a, b, c, d) {
    const n = d.length;
    const cp = new Array(n);
    const dp = new Array(n);
    const x = new Array(n);
    cp[0] = C.div(c, b[0]);
    dp[0] = C.div(d[0], b[0]);
    for (let i = 1; i < n; i++) {
        let den = C.sub(b[i], C.mul(a, cp[i - 1]));
        cp[i] = C.div(c, den);
        dp[i] = C.div(C.sub(d[i], C.mul(a, dp[i - 1])), den);
    }
    x[n - 1] = dp[n - 1];
    for (let i = n - 2; i >= 0; i--) {
        x[i] = C.sub(dp[i], C.mul(cp[i], x[i + 1]));
    }
    return x;
}

function runSimulation() {

    // ==========================================
    //  Inicialização dos Vetores
    // ==========================================
    let uin = new Array(N); // Vetor u interno
    const Vvec = new Array(N);

    for (let i = 0; i < N; i++) {
        let zi = (i + 1) * dz;
        
        // Potencial V 
        let vVal = Math.exp(-Math.pow(zi - 250, 2) / 30000);
        Vvec[i] = { re: vVal, im: 0 };

        // Condição Inicial u0 
        // u0 = exp(1i*z/4) * exp(-(z-250)^2/1500)
        let amp = Math.exp(-Math.pow(zi - 250, 2) / 1500);
        let angle = zi / 4;
        uin[i] = { 
            re: amp * Math.cos(angle), 
            im: amp * Math.sin(angle) 
        };
    }

    // ==========================================
    //  Matrizes Crank-Nicolson (Coeficientes)
    // ==========================================
    // CoefM = (1i * b * dx) / (2 * dz^2)
    const coefM = { re: 0, im: (b * dx) / (2 * dz * dz) };

    // No algoritmo de Thomas para A*uin_new = B*uin_old:
    // d_old = B * uin_old
    // Matriz A é I - M - (i*dx/2)*V
    const mainDiagA = new Array(N);
    const offDiagA = C.scale(coefM, -1); // -1 do termo da segunda derivada

    for (let i = 0; i < N; i++) {
        // Termo de V: (i * dx / 2) * V
        let idv = { re: -Vvec[i].im * (dx / 2), im: Vvec[i].re * (dx / 2) };
        
        // b_ii = 1 - (-2 * coefM) - idv = 1 + 2*coefM - idv
        mainDiagA[i] = C.sub(C.add({ re: 1, im: 0 }, C.scale(coefM, 2)), idv);
    }

    // ==========================================
    //  Evolução em X (Loop de Marcha)
    // ==========================================
    for (let n = 0; n < nx; n++) {
        // 1. Guardar magnitude atual para renderização
        let magProfile = new Float32Array(nz);
        magProfile[0] = 0;
        magProfile[nz - 1] = 0;
        for (let i = 0; i < N; i++) magProfile[i + 1] = C.mag(uin[i]);
        staticSolution.push(magProfile);

        // 2. Calcular Lado Direito: B * uin
        // B = I + M + (i*dx/2)*V
        let rhs = new Array(N);
        for (let i = 0; i < N; i++) {
            let idv = { re: -Vvec[i].im * (dx / 2), im: Vvec[i].re * (dx / 2) };
            
            // Termo central: (1 - 2*coefM + idv) * u[i]
            let centralB = C.add(C.sub({ re: 1, im: 0 }, C.scale(coefM, 2)), idv);
            let term = C.mul(centralB, uin[i]);
            
            // Vizinhos (M): coefM * (u[i+1] + u[i-1])
            if (i > 0) term = C.add(term, C.mul(coefM, uin[i - 1]));
            if (i < N - 1) term = C.add(term, C.mul(coefM, uin[i + 1]));
            
            rhs[i] = term;
        }

        // 3. Resolver A * u_next = rhs usando Thomas
        uin = solveThomas(offDiagA, mainDiagA, offDiagA, rhs);
    }

    document.getElementById('status').innerText = "Simulando propagação...";
    animate();
}

/** * PARTE DO PLOT (CANVAS)
 */
let currentR = 0;
function animate() {
    if (currentR >= staticSolution.length) return;

    const column = staticSolution[currentR];
    const scaleX = canvas.width / nx;
    const scaleY = canvas.height / nz;

    for (let z = 0; z < nz; z += 2) { // z+=2 para performance
        let val = column[z] * 255; 
        
        // Mapeamento de cor similar ao MATLAB 'imagesc'
        ctx.fillStyle = `rgb(${val}, ${val * 0.5}, ${255 - val})`;
        ctx.fillRect(currentR * scaleX, canvas.height - (z * scaleY), scaleX + 1, scaleY * 2);
    }

    currentR++;
    requestAnimationFrame(animate);
}

setTimeout(runSimulation, 100);
