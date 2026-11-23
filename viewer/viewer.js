(async function () {
  'use strict';
  const data = window.JACKET_LOAD_DATA;
  const chart = document.getElementById('chart');
  const slider = document.getElementById('phase');
  const play = document.getElementById('play');
  const select = document.getElementById('viewpoint');
  const error = document.getElementById('error');
  if (!data || !window.Plotly) {
    error.hidden = false;
    error.textContent = 'The interactive graph could not load. Reload the page or download the AVI.';
    return;
  }
  const n = data.headings.length;
  const axes = [...data.headings, 360];
  const frameCount = data.phases.length;
  const surfaces = data.loads.map(rows => {
    const out = rows.map(row => [...row, row[0]]);
    out.push([...out[0]]);
    return out;
  });
  const peaks = data.loads.map(rows => {
    let value = -Infinity, i = 0, j = 0;
    rows.forEach((row, r) => row.forEach((v, c) => {
      if (v > value) { value = v; i = r; j = c; }
    }));
    return {value, wave: data.headings[i], current: data.headings[j]};
  });
  const ceiling = Math.ceil(Math.max(...peaks.map(p => p.value)));
  const cameras = {
    oblique: {eye:{x:1.45,y:1.65,z:1.25},up:{x:0,y:0,z:1}},
    top: {eye:{x:0,y:0,z:2.4},up:{x:0,y:1,z:0}},
    current: {eye:{x:0,y:-2.5,z:0.05},up:{x:0,y:0,z:1}},
    wave: {eye:{x:2.5,y:0,z:0.05},up:{x:0,y:0,z:1}}
  };
  let frame = 0, playing = false, busy = false, pending = null;
  let animationStart = 0;
  const first = peaks[0];
  const traces = [
    {type:'surface', x:axes, y:axes, z:surfaces[0], colorscale:'Viridis',
      cmin:0,cmax:ceiling,colorbar:{title:{text:'Load (MN)'},len:.62,thickness:14},
      hovertemplate:'Current %{x:.0f}°<br>Wave %{y:.0f}°<br>Load %{z:.3f} MN<extra></extra>'},
    {type:'scatter3d',mode:'markers',x:[first.current],y:[first.wave],z:[first.value],
      marker:{size:5,color:'#df3b47',line:{color:'#fff',width:2}},name:'Sampled peak',showlegend:false,
      hovertemplate:'Sampled peak at this phase<br>Current %{x:.0f}°<br>Wave %{y:.0f}°<br>Load %{z:.3f} MN<extra></extra>'}
  ];
  const axis = {range:[0,360],tickvals:[0,90,180,270,360],gridcolor:'#e1e7ef',zerolinecolor:'#bdc9d9',tickfont:{size:12}};
  try {
    await Plotly.newPlot(chart,traces,{
      margin:{l:10,r:25,t:12,b:12},paper_bgcolor:'#fff',font:{family:'system-ui, sans-serif',size:14,color:'#243a52'},
      uirevision:'retain-view',scene:{uirevision:'retain-view',dragmode:'orbit',camera:cameras.oblique,
        xaxis:{...axis,title:{text:'Current heading (°)'}},
        yaxis:{...axis,title:{text:'Wave heading (°)'}},
        zaxis:{range:[0,ceiling],title:{text:'Horizontal load (MN)'},gridcolor:'#e1e7ef',tickfont:{size:12}},
        aspectmode:'manual',aspectratio:{x:1,y:1,z:.7}}
    },{responsive:true,displaylogo:false,scrollZoom:true,displayModeBar:false});
  } catch (err) {
    error.hidden=false;error.textContent='Your browser could not display the 3D graph. Try another browser or download the AVI.';
    return;
  }
  function labels(index) {
    const phase=data.phases[index], peak=peaks[index];
    document.getElementById('phase-value').textContent=phase.toFixed(1)+'°';
    document.getElementById('time-value').textContent=(phase/360*data.period).toFixed(2)+' / '+data.period.toFixed(2)+' s';
    document.getElementById('peak').textContent='Sampled peak at this phase: '+peak.value.toFixed(3)+' MN · wave '+peak.wave+'° · current '+peak.current+'°';
    slider.value=index;slider.setAttribute('aria-valuetext',phase.toFixed(1)+' degrees');
  }
  async function show(index) {
    pending=index;
    if(busy) return;
    busy=true;
    try {
      while(pending !== null) {
        const next=pending;pending=null;
        const peak=peaks[next];
        await Plotly.restyle(chart,{z:[surfaces[next],[peak.value]],x:[axes,[peak.current]],y:[axes,[peak.wave]]},[0,1]);
        frame=next;labels(next);
      }
    } finally {busy=false;}
  }
  function stop() {
    playing=false;play.textContent='Play cycle';play.setAttribute('aria-pressed','false');
  }
  function tick(now) {
    if(!playing) return;
    const next=Math.floor(((now-animationStart)%(data.period*1000))/(data.period*1000)*frameCount);
    if(next !== frame) show(next);
    requestAnimationFrame(tick);
  }
  play.addEventListener('click',()=>{
    if(playing){stop();return;}
    playing=true;animationStart=performance.now()-frame/frameCount*data.period*1000;
    play.textContent='Pause';play.setAttribute('aria-pressed','true');requestAnimationFrame(tick);
  });
  slider.max=frameCount-1;
  slider.addEventListener('input',()=>{stop();show(Number(slider.value));});
  select.addEventListener('change',()=>Plotly.relayout(chart,{'scene.camera':cameras[select.value]}));
  document.addEventListener('visibilitychange',()=>{if(document.hidden)stop();});
  labels(0);
})();
