(function () { 
   // --- Parâmetros Fixos da Grade ---
    const z_plot_min = -800.0;
    const z_plot_max = 200.0;
    const z_grid_max = 600.0;
    const dz = 0.5;
    const r_max = 10000.0;
    const dr = 5.0;

    const Nz = Math.floor((z_grid_max - z_plot_min) / dz) + 1;
    const Nr = Math.floor(r_max / dr) + 1;
    const idx_min = 0;
    const idx_max = Math.floor((z_plot_max - z_plot_min) / dz) + 1;
    const plot_Nz = idx_max - idx_min;

    let z_array = new Float64Array(Nz);
    for(let i=0; i<Nz; i++) z_array[i] = z_plot_min + i * dz;

    // Lógica Visual e Drag
    const bolinha = document.getElementById('source-dot1');
    const canvas = document.getElementById('simCanvas1');
    const ctx = canvas.getContext('2d');
    let arrastando = false;

    // Conversão Coordenada -> Pixel
    function zToY(z_val) {
        return ((z_val - z_plot_min) / (z_plot_max - z_plot_min)) * canvas.height;
    }
    function yToZ(yPx) {
        return (yPx / canvas.height) * (z_plot_max - z_plot_min) + z_plot_min;
    }

    function atualizarPosicaoBolinha() {
        let zs = parseFloat(document.getElementById('s1_zs').value);
        if(isNaN(zs)) return;
        // Restringe fonte à coluna d'água (<= 0) e ao min (-800)
        if (zs < z_plot_min) zs = z_plot_min;
        if (zs > 0) zs = 0;
    
        bolinha.style.top = zToY(zs) + 'px';
    }

    bolinha.addEventListener('mousedown', (e) => {
        arrastando = true; e.preventDefault();
    });
    document.addEventListener('mousemove', (e) => {
        if (!arrastando) return;
        const rect = canvas.getBoundingClientRect();
        let y = e.clientY - rect.top;
    
        let z_lim_max = zToY(0); // Não passa do fundo do mar (0m)
        if (y < 0) y = 0;
        if (y > z_lim_max) y = z_lim_max;
    
        bolinha.style.top = y + 'px';
        document.getElementById('s1_zs').value = yToZ(y).toFixed(1);
    });
    document.addEventListener('mouseup', () => {
        if (arrastando) { arrastando = false; iniciarSimulacao1(); }
    });

    // --- Algoritmo de Thomas Otimizado (Números Complexos em Arrays Primitivos 64bit) ---
    function spsolve_tdma_inline(a_r, a_i, b_r, b_i, c_r, c_i, d_r, d_i) {
        let n = d_r.length;
        let cp_r = new Float64Array(n), cp_i = new Float64Array(n);
        let dp_r = new Float64Array(n), dp_i = new Float64Array(n);
        let x_r = new Float64Array(n), x_i = new Float64Array(n);

        let den = b_r[0]*b_r[0] + b_i[0]*b_i[0];
        cp_r[0] = (c_r[0]*b_r[0] + c_i[0]*b_i[0])/den;
        cp_i[0] = (c_i[0]*b_r[0] - c_r[0]*b_i[0])/den;

        dp_r[0] = (d_r[0]*b_r[0] + d_i[0]*b_i[0])/den;
        dp_i[0] = (d_i[0]*b_r[0] - d_r[0]*b_i[0])/den;

        for (let i = 1; i < n; i++) {
            let axc_r = a_r[i]*cp_r[i-1] - a_i[i]*cp_i[i-1];
            let axc_i = a_r[i]*cp_i[i-1] + a_i[i]*cp_r[i-1];
            let m_r = b_r[i] - axc_r;
            let m_i = b_i[i] - axc_i;

            den = m_r*m_r + m_i*m_i;

            if (i < n - 1) {
                cp_r[i] = (c_r[i]*m_r + c_i[i]*m_i)/den;
                cp_i[i] = (c_i[i]*m_r - c_r[i]*m_i)/den;
            }

            let axd_r = a_r[i]*dp_r[i-1] - a_i[i]*dp_i[i-1];
            let axd_i = a_r[i]*dp_i[i-1] + a_i[i]*dp_r[i-1];
            let num_r = d_r[i] - axd_r;
            let num_i = d_i[i] - axd_i;

            dp_r[i] = (num_r*m_r + num_i*m_i)/den;
            dp_i[i] = (num_i*m_r - num_r*m_i)/den;
        }

        x_r[n-1] = dp_r[n-1];
        x_i[n-1] = dp_i[n-1];

        for (let i = n - 2; i >= 0; i--) {
            let cxd_r = cp_r[i]*x_r[i+1] - cp_i[i]*x_i[i+1];
            let cxd_i = cp_r[i]*x_i[i+1] + cp_i[i]*x_r[i+1];
            x_r[i] = dp_r[i] - cxd_r;
            x_i[i] = dp_i[i] - cxd_i;
        }
        return {r: x_r, i: x_i};
    }

    // --- Renderização Eficiente (Escala de Magma Reversa) ---
    function getMagmaRColor(val) {
        if (val <= 34) return [255, 255, 255];
        if (val >= 48) return [0, 0, 0];
    
        const t = (val - 34) / (48 - 34);
        const stops = [
            {t: 0.0, c: [255, 255, 255]},
            {t: 0.2, c: [252, 253, 191]},
            {t: 0.4, c: [252, 137, 97]},
            {t: 0.6, c: [183, 55, 121]},
            {t: 0.8, c: [81, 18, 124]},
            {t: 1.0, c: [0, 0, 0]}
        ];
    
        for (let i = 0; i < stops.length - 1; i++) {
            if (t >= stops[i].t && t <= stops[i+1].t) {
                const range = stops[i+1].t - stops[i].t;
                const lt = (t - stops[i].t) / range;
                const r = Math.round(stops[i].c[0] + lt * (stops[i+1].c[0] - stops[i].c[0]));
                const g = Math.round(stops[i].c[1] + lt * (stops[i+1].c[1] - stops[i].c[1]));
                const b = Math.round(stops[i].c[2] + lt * (stops[i+1].c[2] - stops[i].c[2]));
                return [r, g, b];
            }
        }
        return [0,0,0];
    }

    // Variáveis Globais da Simulação
    let imgDataBuffer, offscreenCanvas, offscreenCtx;
    let current_ir = 0;
    let U_final;
    let psi_r, psi_i;
    let A_sub_r, A_sub_i, A_diag_r, A_diag_i, A_sup_r, A_sup_i;
    let B_sub_r, B_sub_i, B_diag_r, B_diag_i, B_sup_r, B_sup_i;
    let runFlag = false;

    // --- Preparação Matemática ---
    function iniciarSimulacao1() {
        runFlag = true;
        document.getElementById('btn-simular').disabled = true;
        document.getElementById('progress-bar').style.display = "block";
        document.getElementById('status-text').innerText = "Montando matrizes...";
    
        const f = parseFloat(document.getElementById('s1_f').value);
        const zs = parseFloat(document.getElementById('s1_zs').value);
        const theta0 = parseFloat(document.getElementById('s1_theta').value);
        const W = parseFloat(document.getElementById('s1_W').value);
        const c_sed = parseFloat(document.getElementById('s1_c_sed').value);
        const rho_sed = parseFloat(document.getElementById('s1_rho_sed').value);
    
        const c0 = 1500.0;
        const k0 = 2 * Math.PI * f / c0;

        let c = new Float64Array(Nz);
        let rho = new Float64Array(Nz);
        let alpha = new Float64Array(Nz);
    
        // Perfil de Meio (Velocidade, Densidade e Atenuação)
        for (let i = 0; i < Nz; i++) {
            let zi = z_array[i];
            if (zi <= 0) {
                c[i] = 1500.0 + 0.016 * (zi + 400.0);
                rho[i] = 1000.0;
            } else {
                c[i] = c_sed;
                rho[i] = rho_sed;
                if (zi > 200) {
                    alpha[i] = 0.5 * Math.pow((zi - 200) / 400.0, 2);
                }
            }
        }

        // Vetores da Matriz IFD
        let rho_mid = new Float64Array(Nz - 1);
        let V_up = new Float64Array(Nz - 1);
        let V_down = new Float64Array(Nz - 1);
        for (let i = 0; i < Nz - 1; i++) {
            rho_mid[i] = (rho[i] + rho[i+1]) / 2.0;
            V_up[i] = rho[i] / rho_mid[i] / (dz * dz);
            V_down[i] = rho[i+1] / rho_mid[i] / (dz * dz);
        }

        let diag_V = new Float64Array(Nz);
        for (let i = 0; i < Nz - 1; i++) {
            diag_V[i] -= V_up[i];
            diag_V[i+1] -= V_down[i];
        }

        // u_diag
        let u_diag_r = new Float64Array(Nz);
        let u_diag_i = new Float64Array(Nz);
        for (let i = 0; i < Nz; i++) {
            let nr = c0 / c[i];
            let ni = alpha[i];
            let n2r = nr*nr - ni*ni;
            let n2i = 2*nr*ni;
            u_diag_r[i] = k0*k0*(n2r - 1) + diag_V[i];
            u_diag_i[i] = k0*k0*n2i;
        }

        // gamma = 1j * dr / (4 * k0)  --> gamma_r = 0, gamma_i = dr / (4*k0)
        let g_i = dr / (4 * k0);

        // Arrays para Solver TDMA
        A_sub_r = new Float64Array(Nz); A_sub_i = new Float64Array(Nz);
        A_diag_r = new Float64Array(Nz); A_diag_i = new Float64Array(Nz);
        A_sup_r = new Float64Array(Nz); A_sup_i = new Float64Array(Nz);

        B_sub_r = new Float64Array(Nz); B_sub_i = new Float64Array(Nz);
        B_diag_r = new Float64Array(Nz); B_diag_i = new Float64Array(Nz);
        B_sup_r = new Float64Array(Nz); B_sup_i = new Float64Array(Nz);

        for (let i = 0; i < Nz; i++) {
            A_diag_r[i] = 1.0 - (-g_i * u_diag_i[i]); A_diag_i[i] = 0.0 - (g_i * u_diag_r[i]);
            B_diag_r[i] = 1.0 + (-g_i * u_diag_i[i]); B_diag_i[i] = 0.0 + (g_i * u_diag_r[i]);

            if (i > 0) {
                A_sub_r[i] = 0.0; A_sub_i[i] = g_i * (-V_down[i-1]);
                B_sub_r[i] = 0.0; B_sub_i[i] = g_i * (V_down[i-1]);
            }
            if (i < Nz - 1) {
                A_sup_r[i] = 0.0; A_sup_i[i] = g_i * (-V_up[i]);
                B_sup_r[i] = 0.0; B_sup_i[i] = g_i * (V_up[i]);
            }
        }

        // Feixe Inicial
        let tilt_rad = theta0 * Math.PI / 180.0;
        psi_r = new Float64Array(Nz);
        psi_i = new Float64Array(Nz);
        for (let i = 0; i < Nz; i++) {
            let z_val = z_array[i];
            let gaus = Math.exp(-Math.pow(z_val - zs, 2) / Math.pow(W, 2));
            let phase = k0 * (z_val - zs) * Math.sin(tilt_rad);
            psi_r[i] = gaus * Math.cos(phase);
            psi_i[i] = gaus * Math.sin(phase);
        }

        // Limpar Canvas Offscreen
        U_final = new Array(Nr);
        if (!offscreenCanvas) {
            offscreenCanvas = document.createElement('canvas');
            offscreenCanvas.width = Nr;
            offscreenCanvas.height = plot_Nz;
            offscreenCtx = offscreenCanvas.getContext('2d');
            imgDataBuffer = offscreenCtx.createImageData(Nr, plot_Nz);
        }
    
        current_ir = 0;
        setTimeout(processarChunk, 10);
    }

    // --- Laço de Marcha Assíncrono para Não Travar o Navegador ---
    function processarChunk() {
        if(!runFlag) return;
        let passos = 50; // Calcula X colunas por frame de renderização
        let end_ir = Math.min(current_ir + passos, Nr);
    
        let dBuffer = imgDataBuffer.data;

        for (let ir = current_ir; ir < end_ir; ir++) {
            // Salva resultado visual
            let dist = ir * dr;
            let p_div = Math.sqrt(dist + 1e-6);
            U_final[ir] = new Float32Array(plot_Nz);
        
            for (let i = 0; i < plot_Nz; i++) {
                let true_idx = i + idx_min;
                let abs_val = Math.sqrt(psi_r[true_idx]*psi_r[true_idx] + psi_i[true_idx]*psi_i[true_idx]);
                let pressao = abs_val / p_div;
                let tl = -20 * Math.log10(pressao + 1e-12);
                U_final[ir][i] = tl;
            
                // Joga direto no buffer de pixels (alta velocidade)
                let color = getMagmaRColor(tl);
                let px = (i * Nr + ir) * 4;
                dBuffer[px] = color[0]; dBuffer[px+1] = color[1]; dBuffer[px+2] = color[2]; dBuffer[px+3] = 255;
            }

            // Multiplicação do RHS (B * psi)
            let rhs_r = new Float64Array(Nz); let rhs_i = new Float64Array(Nz);
            for(let i = 0; i < Nz; i++){
                let tr = B_diag_r[i]*psi_r[i] - B_diag_i[i]*psi_i[i];
                let ti = B_diag_r[i]*psi_i[i] + B_diag_i[i]*psi_r[i];
                if (i > 0) {
                    tr += B_sub_r[i]*psi_r[i-1] - B_sub_i[i]*psi_i[i-1];
                    ti += B_sub_r[i]*psi_i[i-1] + B_sub_i[i]*psi_r[i-1];
                }
                if (i < Nz - 1) {
                    tr += B_sup_r[i]*psi_r[i+1] - B_sup_i[i]*psi_i[i+1];
                    ti += B_sup_r[i]*psi_i[i+1] + B_sup_i[i]*psi_r[i+1];
                }
                rhs_r[i] = tr; rhs_i[i] = ti;
            }

            // Resolve Sistema (A \ RHS)
            let sol = spsolve_tdma_inline(A_sub_r, A_sub_i, A_diag_r, A_diag_i, A_sup_r, A_sup_i, rhs_r, rhs_i);
            psi_r = sol.r; psi_i = sol.i;
        
            // Condição de Fronteira Superior (Superfície livre)
            psi_r[0] = 0.0; psi_i[0] = 0.0;
        }

        current_ir = end_ir;
    
        // Atualiza Progresso e Canvas Principal
        let pct = (current_ir / Nr) * 100;
        document.getElementById('progress-fill').style.width = pct + "%";
        document.getElementById('status-text').innerText = `Processando Range: ${(current_ir * dr / 1000).toFixed(1)} km`;

        // Desenha na tela (Esticando a imagem offscreen para caber no 800x400)
        offscreenCtx.putImageData(imgDataBuffer, 0, 0);
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.drawImage(offscreenCanvas, 0, 0, canvas.width, canvas.height);
    
        // Desenha Linha do Fundo do Mar (z=0)
        let linhaFundoY = zToY(0);
        ctx.beginPath();
        ctx.setLineDash([5, 5]);
        ctx.moveTo(0, linhaFundoY);
        ctx.lineTo(canvas.width, linhaFundoY);
        ctx.strokeStyle = "black";
        ctx.lineWidth = 1;
        ctx.stroke();
        ctx.setLineDash([]);

        if (current_ir < Nr) {
            requestAnimationFrame(processarChunk);
        } else {
            document.getElementById('status-text').innerText = "Simulação concluída! (10 km atingidos)";
            document.getElementById('btn-simular').disabled = false;
            runFlag = false;
        }
    }

    // Inicialização
    window.onload = () => {
        atualizarPosicaoBolinha();
        iniciarSimulacao1();
    };
})();