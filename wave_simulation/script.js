const transformations1 = document.querySelector('#transformations1');
const transformations2 = document.querySelector('#transformations2');
transformations1.addEventListener('load', () => {
    const svgDoc = transformations1.contentDocument
    const translation_bar = transformations1.contentDocument.querySelector('#Translation_bar_slider');
    const scale_bar = transformations1.contentDocument.querySelector('#Scale_bar_slider');
    const rotation_bar = transformations1.contentDocument.querySelector('#Rotation_bar_slider');
    const skew_bar = transformations1.contentDocument.querySelector('#Skew_bar_slider');
    const rotation_button = transformations1.contentDocument.querySelector('#Rotation_button_press');
    const skew_button = transformations1.contentDocument.querySelector('#Skew_button_press');
    const figure = transformations1.contentDocument.querySelector('#figure');

    let pressed = false;
    let dx = 0;
    let dy = 0;
    let y_pos = 0;
    let lastY = 0;
    let x_pos = 0;
    let lastX = 0;
    let bar = 0;
    let angle_rotation = 0;
    let angle_skew = 0;


    // Posição de figure
    let x1 = 300;
    let x2 = x1+40;
    let x3 = x1-40;
    let y1 = 180;
    let y2 = y1+150;
    let y3 = y1+50;
    let scale = 1;

    translation_bar.setAttribute("x",x1);
    scale_bar.setAttribute("x",x1);
    
    // Controle da barra de Transladar
    translation_bar.addEventListener('mousedown',(e)=>{
        pressed = true;
        bar = 0;
        lastX = e.offsetX;
    });

    function checkboundary_translation(){
        if(parseFloat(translation_bar.getAttribute("x") || 0) > 700){
            translation_bar.setAttribute("x",700);
        }
        if(parseFloat(translation_bar.getAttribute("x") || 0) < -120){
            translation_bar.setAttribute("x",-120);
        }
    }

    // Controle da barra de Escalar
    scale_bar.addEventListener('mousedown',(e)=>{
        pressed = true;
        bar = 1;
        lastX = e.offsetX;
    });
 
    function checkboundary_scale(){
        if(parseFloat(scale_bar.getAttribute("x") || 0) > 700){
            scale_bar.setAttribute("x",700);
        }
        if(parseFloat(scale_bar.getAttribute("x") || 0) < -120){
            scale_bar.setAttribute("x",-120);
        }
    }
    
    // Controle da barra de Rotação
    rotation_bar.addEventListener('mousedown',(e)=>{
        pressed = true;
        bar = 2;
        lastY = e.offsetY;
    });

    function checkboundary_rotation(){
        if(parseFloat(rotation_bar.getAttribute("y") || 0) < 75){
            rotation_bar.setAttribute("y",75);
        }
        if(parseFloat(rotation_bar.getAttribute("y") || 0) > 480){
            rotation_bar.setAttribute("y",480);
        }
    }


    // Controle da barra de Skew
    skew_bar.addEventListener('mousedown',(e)=>{
        pressed = true;
        bar = 3;
        lastY = e.offsetY;
    });

    function checkboundary_skew(){
        if(parseFloat(skew_bar.getAttribute("y") || 0) < 75){
            skew_bar.setAttribute("y",75);
        }
        if(parseFloat(skew_bar.getAttribute("y") || 0) > 480){
            skew_bar.setAttribute("y",480);
        }
    }


    // Mudança para todos os conotrles
    svgDoc.addEventListener('mousemove',(e)=>{
        if (!pressed) return;
        e.preventDefault();
        
        dx = e.offsetX - lastX;
        lastX = e.offsetX;
        
        dy = e.offsetY - lastY;
        lastY = e.offsetY;
        
        if(bar==0){
        x_pos = parseFloat(translation_bar.getAttribute("x") || 0);
        translation_bar.setAttribute("x",x_pos+dx);
        x1 = parseFloat(translation_bar.getAttribute("x"));
        checkboundary_translation();
        figure_render(x1,y1,scale);
        }
        else if(bar==1){
        x_pos = parseFloat(scale_bar.getAttribute("x") || 0);    
        scale_bar.setAttribute("x",x_pos+dx);
        scale = parseFloat(scale_bar.getAttribute("x"))/600+1/2;
        checkboundary_scale();   
        figure_render(x1,y1,scale);         
        }
        else if(bar==2){
        y_pos = parseFloat(rotation_bar.getAttribute("y") || 0);    
        rotation_bar.setAttribute("y",y_pos+dy);
        checkboundary_rotation(); 
        let y_barra = parseFloat(rotation_bar.getAttribute("y"));

        let angle_rotation = (y_barra - 75) * (360 / (480 - 75));

        figure.setAttribute("transform", `translate(${x1}, ${y2}) skewX(${angle_skew}) translate(${-x1}, ${-y2}) rotate(${angle_rotation}, ${x1}, ${y1})`);
        }
        else if(bar==3){
        y_pos = parseFloat(skew_bar.getAttribute("y") || 0);    
        skew_bar.setAttribute("y",y_pos+dy);
        checkboundary_skew();    
        
        let y_val = parseFloat(skew_bar.getAttribute("y"));
        angle_skew = (y_val - 75) * (90 / (480 - 75)) - 45;
        figure.setAttribute("transform", `translate(${x1}, ${y2}) skewX(${angle_skew}) translate(${-x1}, ${-y2}) rotate(${angle_rotation}, ${x1}, ${y1})`);
        }
    });

    svgDoc.addEventListener('mouseup',()=>{
        pressed = false;
    });

    // Atualizar figure
    function figure_render(x1,y1,scale){
        x2 = x1+40*scale;
        x3 = x1-40*scale;
        y2 = y1+150*scale;
        y3 = y1+50*scale;        
        figure.setAttribute('points', `${x1},${y1} ${x2},${y1} ${x2},${y2} ${x1},${y2} ${x1},${y3} ${x3},${y3}`)
    }
    }); 
transformations2.addEventListener('load', () => {
    const svgDoc2 = transformations2.contentDocument;
    const canvas = svgDoc2.querySelector('svg');

    let pontosControle = [];
    let pontoSelecionado = null; // Guarda o índice do ponto que está sendo arrastado

    const cubicaPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
    cubicaPath.setAttribute("fill", "none");
    cubicaPath.setAttribute("stroke", "yellow");
    cubicaPath.setAttribute("stroke-width", "3");
    canvas.appendChild(cubicaPath);

    // Função para atualizar a string do Path
    function atualizarBezier() {
        if (pontosControle.length === 4) {
            const [p1, p2, p3, p4] = pontosControle;
            const d = `M ${p1.x} ${p1.y} C ${p2.x} ${p2.y}, ${p3.x} ${p3.y}, ${p4.x} ${p4.y}`;
            cubicaPath.setAttribute("d", d);
        }
    }

    canvas.addEventListener('mousedown', (e) => {
        if (e.target.tagName === 'circle') {
            pontoSelecionado = parseInt(e.target.getAttribute('data-index'));
            return;
        }

        if (pontosControle.length < 4) {
            let ponto = canvas.createSVGPoint();
            ponto.x = e.clientX;
            ponto.y = e.clientY;
            const coords = ponto.matrixTransform(canvas.getScreenCTM().inverse());

            const idx = pontosControle.length;
            pontosControle.push({ x: coords.x, y: coords.y });

            const new_circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
            new_circle.setAttribute("cx", coords.x);
            new_circle.setAttribute("cy", coords.y);
            new_circle.setAttribute("r", 10);
            new_circle.setAttribute("fill", "white");
            new_circle.setAttribute("stroke", "black");
            new_circle.setAttribute("data-index", idx); 
            new_circle.style.cursor = "move";
            
            canvas.appendChild(new_circle);
            atualizarBezier();
        }
    });

    canvas.addEventListener('mousemove', (e) => {
        if (pontoSelecionado === null) return;

        let ponto = canvas.createSVGPoint();
        ponto.x = e.clientX;
        ponto.y = e.clientY;
        const coords = ponto.matrixTransform(canvas.getScreenCTM().inverse());

        // Atualiza a posição no Array
        pontosControle[pontoSelecionado].x = coords.x;
        pontosControle[pontoSelecionado].y = coords.y;

        // Atualiza a posição visual do Círculo
        const circulo = canvas.querySelector(`circle[data-index="${pontoSelecionado}"]`);
        circulo.setAttribute("cx", coords.x);
        circulo.setAttribute("cy", coords.y);

        atualizarBezier();
    });

    canvas.addEventListener('mouseup', () => {
        pontoSelecionado = null;
    });
});