/* ============================================
   iimjobs Experiment Hub — Application Logic
   ============================================ */

// ========================
// NAVIGATION
// ========================
document.addEventListener('DOMContentLoaded', () => {
    initNavigation();
    initThemeToggle();
    initMobileMenu();
    initChecklistProgress();
    initSQLBlocks();
    initRegistryData();
    initDashboardCharts();
});

function initNavigation() {
    const navItems = document.querySelectorAll('.nav-item');
    navItems.forEach(item => {
        item.addEventListener('click', (e) => {
            e.preventDefault();
            const sectionId = item.getAttribute('data-section');
            switchSection(sectionId);
            // Close mobile sidebar
            document.getElementById('sidebar').classList.remove('open');
        });
    });
    // Handle hash on load
    const hash = window.location.hash.replace('#', '');
    if (hash) switchSection(hash);
}

function switchSection(sectionId) {
    document.querySelectorAll('.content-section').forEach(s => s.classList.remove('active'));
    document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));

    const section = document.getElementById(sectionId);
    const navItem = document.querySelector(`[data-section="${sectionId}"]`);
    if (section) section.classList.add('active');
    if (navItem) navItem.classList.add('active');

    window.location.hash = sectionId;
    window.scrollTo({ top: 0 });

    // Re-render charts if switching to dashboard
    if (sectionId === 'dashboard') {
        setTimeout(initDashboardCharts, 100);
    }
}

// ========================
// THEME TOGGLE
// ========================
function initThemeToggle() {
    const stored = localStorage.getItem('theme') || 'dark';
    document.documentElement.setAttribute('data-theme', stored);

    document.getElementById('themeToggle').addEventListener('click', toggleTheme);
    const mobileToggle = document.getElementById('themeToggleMobile');
    if (mobileToggle) mobileToggle.addEventListener('click', toggleTheme);
}

function toggleTheme() {
    const current = document.documentElement.getAttribute('data-theme');
    const next = current === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    localStorage.setItem('theme', next);
    const mobileToggle = document.getElementById('themeToggleMobile');
    if (mobileToggle) mobileToggle.textContent = next === 'dark' ? '🌙' : '☀️';
    // Redraw charts with new theme
    if (document.getElementById('dashboard').classList.contains('active')) {
        setTimeout(initDashboardCharts, 50);
    }
}

// ========================
// MOBILE MENU
// ========================
function initMobileMenu() {
    document.getElementById('hamburger').addEventListener('click', () => {
        document.getElementById('sidebar').classList.toggle('open');
    });
    // Close on outside click
    document.getElementById('mainContent').addEventListener('click', () => {
        document.getElementById('sidebar').classList.remove('open');
    });
}

// ========================
// STATISTICAL FUNCTIONS
// ========================

// Standard normal CDF (Abramowitz & Stegun approximation)
function normalCDF(x) {
    const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x) / Math.sqrt(2);
    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return 0.5 * (1.0 + sign * y);
}

// Inverse normal CDF (rational approximation)
function normalInv(p) {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    if (p === 0.5) return 0;

    const a = [
        -3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285104469687e+02, 1.383577518672690e+02,
        -3.066479806614716e+01, 2.506628277459239e+00
    ];
    const b = [
        -5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01,
        -1.328068155288572e+01
    ];
    const c = [
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
        4.374664141464968e+00, 2.938163982698783e+00
    ];
    const d = [
        7.784695709041462e-03, 3.224671290700398e-01,
        2.445134137142996e+00, 3.754408661907416e+00
    ];

    const p_low = 0.02425;
    const p_high = 1 - p_low;
    let q, r;

    if (p < p_low) {
        q = Math.sqrt(-2 * Math.log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    } else if (p <= p_high) {
        q = p - 0.5;
        r = q * q;
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    } else {
        q = Math.sqrt(-2 * Math.log(1 - p));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
}

// Chi-square CDF (1 degree of freedom)
function chiSquareCDF1(x) {
    if (x <= 0) return 0;
    return 2 * normalCDF(Math.sqrt(x)) - 1;
}

// ========================
// SAMPLE SIZE CALCULATOR
// ========================
function calculateSampleSize() {
    const baseline = parseFloat(document.getElementById('ssBaseline').value) / 100;
    const mde = parseFloat(document.getElementById('ssMDE').value) / 100;
    const confidence = parseFloat(document.getElementById('ssConfidence').value) / 100;
    const power = parseFloat(document.getElementById('ssPower').value) / 100;
    const dailyTraffic = parseInt(document.getElementById('ssDailyTraffic').value);

    if (isNaN(baseline) || isNaN(mde) || isNaN(confidence) || isNaN(power) || isNaN(dailyTraffic)) {
        alert('Please fill in all fields with valid numbers.');
        return;
    }

    const p1 = baseline;
    const p2 = baseline * (1 + mde);

    const zAlpha = normalInv(1 - (1 - confidence) / 2);
    const zBeta = normalInv(power);

    // Two-proportion z-test formula
    const n = Math.ceil(
        Math.pow(zAlpha + zBeta, 2) * (p1 * (1 - p1) + p2 * (1 - p2)) /
        Math.pow(p2 - p1, 2)
    );

    const totalN = n * 2;
    const days = Math.ceil(n / dailyTraffic);

    document.getElementById('ssResultN').textContent = n.toLocaleString();
    document.getElementById('ssResultTotal').textContent = totalN.toLocaleString();
    document.getElementById('ssResultDuration').textContent = `${days} days`;
    document.getElementById('ssResultTreatment').textContent = `${(p2 * 100).toFixed(2)}%`;
    document.getElementById('ssFormula').style.display = 'block';
}

// ========================
// SRM CHECKER
// ========================
function checkSRM() {
    const control = parseInt(document.getElementById('srmControl').value);
    const treatment = parseInt(document.getElementById('srmTreatment').value);
    const expectedPct = parseFloat(document.getElementById('srmExpected').value) / 100;

    if (isNaN(control) || isNaN(treatment) || isNaN(expectedPct)) {
        alert('Please fill in all fields.');
        return;
    }

    const total = control + treatment;
    const expectedControl = total * expectedPct;
    const expectedTreatment = total * (1 - expectedPct);

    const chiSquare = Math.pow(control - expectedControl, 2) / expectedControl +
                      Math.pow(treatment - expectedTreatment, 2) / expectedTreatment;

    const pValue = 1 - chiSquareCDF1(chiSquare);

    document.getElementById('srmChi').textContent = chiSquare.toFixed(4);
    document.getElementById('srmPval').textContent = pValue < 0.0001 ? pValue.toExponential(4) : pValue.toFixed(6);

    const verdictEl = document.getElementById('srmVerdict');
    const verdictText = document.getElementById('srmVerdictText');

    if (pValue < 0.01) {
        verdictText.textContent = '🚨 SRM DETECTED — Randomization Invalid';
        verdictText.style.color = 'var(--danger)';
        verdictEl.style.borderColor = 'var(--danger)';
        verdictEl.style.background = 'var(--danger-bg)';
    } else {
        verdictText.textContent = '✅ Randomization Valid — No SRM Detected';
        verdictText.style.color = 'var(--success)';
        verdictEl.style.borderColor = 'var(--success)';
        verdictEl.style.background = 'var(--success-bg)';
    }
}

// ========================
// LIFT CALCULATOR
// ========================
function calculateLift() {
    const controlConv = parseInt(document.getElementById('liftControlConv').value);
    const controlN = parseInt(document.getElementById('liftControlN').value);
    const treatConv = parseInt(document.getElementById('liftTreatConv').value);
    const treatN = parseInt(document.getElementById('liftTreatN').value);

    if (isNaN(controlConv) || isNaN(controlN) || isNaN(treatConv) || isNaN(treatN)) {
        alert('Please fill in all fields.');
        return;
    }

    const pC = controlConv / controlN;
    const pT = treatConv / treatN;
    const absLift = pT - pC;
    const relLift = pC === 0 ? 0 : (pT - pC) / pC;

    // Pooled SE for z-test
    const pPooled = (controlConv + treatConv) / (controlN + treatN);
    const se = Math.sqrt(pPooled * (1 - pPooled) * (1/controlN + 1/treatN));
    const z = se === 0 ? 0 : absLift / se;
    const pValue = 2 * (1 - normalCDF(Math.abs(z)));

    // 95% CI on absolute lift
    const seUnpooled = Math.sqrt(pC*(1-pC)/controlN + pT*(1-pT)/treatN);
    const ciLow = absLift - 1.96 * seUnpooled;
    const ciHigh = absLift + 1.96 * seUnpooled;

    document.getElementById('liftCR').textContent = `${(pC * 100).toFixed(3)}%`;
    document.getElementById('liftTR').textContent = `${(pT * 100).toFixed(3)}%`;
    document.getElementById('liftAbs').textContent = `${(absLift * 100).toFixed(3)}pp`;
    document.getElementById('liftRel').textContent = `${(relLift * 100).toFixed(2)}%`;
    document.getElementById('liftZ').textContent = z.toFixed(4);
    document.getElementById('liftPval').textContent = pValue < 0.0001 ? pValue.toExponential(4) : pValue.toFixed(6);
    document.getElementById('liftCI').textContent = `[${(ciLow * 100).toFixed(3)}pp, ${(ciHigh * 100).toFixed(3)}pp]`;

    const verdictEl = document.getElementById('liftVerdict');
    const verdictText = document.getElementById('liftVerdictText');

    if (pValue < 0.05) {
        verdictText.textContent = `✅ Statistically Significant (p < 0.05)`;
        verdictText.style.color = 'var(--success)';
        verdictEl.style.borderColor = 'var(--success)';
        verdictEl.style.background = 'var(--success-bg)';
    } else {
        verdictText.textContent = `⬜ Not Significant (p ≥ 0.05)`;
        verdictText.style.color = 'var(--warning)';
        verdictEl.style.borderColor = 'var(--warning)';
        verdictEl.style.background = 'var(--warning-bg)';
    }
}

// ========================
// SQL COLLAPSIBLE BLOCKS
// ========================
function initSQLBlocks() {
    // Open first block by default
    const firstBlock = document.querySelector('.sql-block');
    if (firstBlock) firstBlock.classList.add('open');
}

function toggleSQL(header) {
    header.parentElement.classList.toggle('open');
}

function copySQL(btn) {
    const code = btn.parentElement.querySelector('code').textContent;
    navigator.clipboard.writeText(code).then(() => {
        const original = btn.textContent;
        btn.textContent = 'Copied!';
        btn.style.background = 'var(--success)';
        setTimeout(() => {
            btn.textContent = original;
            btn.style.background = '';
        }, 1500);
    });
}

// ========================
// CHECKLIST PROGRESS
// ========================
function initChecklistProgress() {
    const checkboxes = document.querySelectorAll('#checklist input[type="checkbox"]');
    const total = checkboxes.length;
    checkboxes.forEach(cb => {
        cb.addEventListener('change', () => {
            const checked = document.querySelectorAll('#checklist input[type="checkbox"]:checked').length;
            const pct = (checked / total) * 100;
            document.getElementById('checklistProgress').style.width = pct + '%';
            document.getElementById('checklistText').textContent = `${checked} / ${total} completed`;
        });
    });
}

// ========================
// EXPERIMENT REGISTRY
// ========================
const defaultExperiments = [
    {
        id: 'EXP-2026-042', hypothesis: 'Simplified apply flow increases applications',
        owner: 'Priya S.', start: '2026-01-15', end: '2026-01-29',
        status: 'Live', metric: 'Application Rate', result: '+9.8%', decision: 'Pending'
    },
    {
        id: 'EXP-2026-039', hypothesis: 'Personalized job recs improve CTR',
        owner: 'Arun K.', start: '2025-12-10', end: '2025-12-24',
        status: 'Completed', metric: 'Job Card CTR', result: '+4.2%', decision: 'Ship'
    },
    {
        id: 'EXP-2026-037', hypothesis: 'Recruiter dashboard redesign increases post rate',
        owner: 'Meera J.', start: '2025-11-20', end: '2025-12-04',
        status: 'Completed', metric: 'Jobs Posted / Recruiter', result: '-1.1%', decision: 'Kill'
    },
    {
        id: 'EXP-2026-035', hypothesis: 'Push notifications boost return visits',
        owner: 'Vikram D.', start: '2025-11-01', end: '2025-11-15',
        status: 'Completed', metric: 'D7 Retention', result: '+2.5%', decision: 'Ship'
    },
    {
        id: 'EXP-2026-044', hypothesis: 'AI-generated cover letters increase apply rate',
        owner: 'Priya S.', start: '2026-02-01', end: '—',
        status: 'Design', metric: 'Application Rate', result: '—', decision: '—'
    },
];

function initRegistryData() {
    const body = document.getElementById('registryBody');
    body.innerHTML = '';
    defaultExperiments.forEach(exp => addRegistryRowData(exp));
}

function addRegistryRowData(data) {
    const body = document.getElementById('registryBody');
    const row = document.createElement('tr');

    const statusOptions = ['Design', 'QA', 'Live', 'Completed', 'Killed'];
    const decisionOptions = ['—', 'Pending', 'Ship', 'Iterate', 'Kill'];

    row.innerHTML = `
        <td><input type="text" value="${data.id}" style="font-family:var(--font-mono);font-weight:600;width:120px"></td>
        <td><input type="text" value="${data.hypothesis}" style="min-width:200px"></td>
        <td><input type="text" value="${data.owner}" style="width:100px"></td>
        <td><input type="date" value="${data.start}"></td>
        <td><input type="${data.end === '—' ? 'text' : 'date'}" value="${data.end}" style="width:110px"></td>
        <td><select>${statusOptions.map(s => `<option${s===data.status?' selected':''}>${s}</option>`).join('')}</select></td>
        <td><input type="text" value="${data.metric}" style="min-width:130px"></td>
        <td><input type="text" value="${data.result}" style="width:80px;font-family:var(--font-mono)"></td>
        <td><select>${decisionOptions.map(d => `<option${d===data.decision?' selected':''}>${d}</option>`).join('')}</select></td>
        <td><button class="btn-delete" onclick="this.closest('tr').remove()">✕</button></td>
    `;
    body.appendChild(row);
}

function addRegistryRow() {
    const newExp = {
        id: `EXP-2026-${String(Math.floor(Math.random()*900)+100).padStart(3,'0')}`,
        hypothesis: '', owner: '', start: new Date().toISOString().slice(0,10),
        end: '—', status: 'Design', metric: '', result: '—', decision: '—'
    };
    addRegistryRowData(newExp);
    // Scroll to bottom of table
    const wrapper = document.querySelector('#registry .table-wrapper');
    wrapper.scrollTop = wrapper.scrollHeight;
}

// ========================
// EXPORT CSV
// ========================
function exportCSV() {
    const rows = document.querySelectorAll('#registryTable tbody tr');
    const headers = ['Experiment ID','Hypothesis','Owner','Start','End','Status','Primary Metric','Result','Decision'];
    let csv = headers.join(',') + '\n';

    rows.forEach(row => {
        const cells = [];
        row.querySelectorAll('input, select').forEach(el => {
            let val = el.value.replace(/"/g, '""');
            cells.push(`"${val}"`);
        });
        // Remove last cell (delete button col)
        csv += cells.join(',') + '\n';
    });

    downloadFile(csv, 'experiment_registry.csv', 'text/csv');
}

// ========================
// EXPORT PDF (Client-side)
// ========================
function exportPDF() {
    const rows = document.querySelectorAll('#registryTable tbody tr');
    const data = [];
    rows.forEach(row => {
        const cells = [];
        row.querySelectorAll('input, select').forEach(el => cells.push(el.value));
        data.push(cells);
    });

    // Build a printable HTML
    let html = `<!DOCTYPE html><html><head><title>Experiment Registry Report</title>
    <style>
        body{font-family:Arial,sans-serif;padding:40px;color:#1a1d2e}
        h1{font-size:20px;margin-bottom:4px}
        p{color:#666;margin-bottom:20px;font-size:12px}
        table{width:100%;border-collapse:collapse;font-size:11px}
        th{background:#1f2937;color:#fff;padding:8px 10px;text-align:left}
        td{padding:7px 10px;border-bottom:1px solid #e0e3ed}
        tr:nth-child(even){background:#f5f6fa}
        .footer{margin-top:30px;font-size:10px;color:#999;text-align:center}
    </style></head><body>
    <h1>iimjobs Experiment Registry</h1>
    <p>Generated: ${new Date().toLocaleString()} | Total experiments: ${data.length}</p>
    <table><thead><tr>
        <th>ID</th><th>Hypothesis</th><th>Owner</th><th>Start</th><th>End</th><th>Status</th><th>Metric</th><th>Result</th><th>Decision</th>
    </tr></thead><tbody>`;

    data.forEach(row => {
        html += '<tr>' + row.map(c => `<td>${escapeHtml(c)}</td>`).join('') + '</tr>';
    });

    html += `</tbody></table>
    <div class="footer">iimjobs Experiment Hub — Confidential</div>
    </body></html>`;

    const printWindow = window.open('', '_blank');
    printWindow.document.write(html);
    printWindow.document.close();
    printWindow.focus();
    setTimeout(() => printWindow.print(), 300);
}

function escapeHtml(str) {
    const div = document.createElement('div');
    div.textContent = str;
    return div.innerHTML;
}

function downloadFile(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}

// ========================
// DASHBOARD CHARTS
// ========================
let chartInstances = {};

function initDashboardCharts() {
    const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
    const gridColor = isDark ? 'rgba(42,45,62,0.6)' : 'rgba(224,227,237,0.8)';
    const textColor = isDark ? '#8b8fa8' : '#5c6080';
    const controlColor = '#4f7cff';
    const treatColor = '#a78bfa';

    const days = ['Day 1','Day 2','Day 3','Day 4','Day 5','Day 6','Day 7','Day 8','Day 9'];

    const defaultOpts = {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
            legend: {
                display: true,
                position: 'top',
                labels: { color: textColor, font: { size: 11, family: 'DM Sans' }, boxWidth: 12, padding: 12 }
            }
        },
        scales: {
            x: { ticks: { color: textColor, font: { size: 10 } }, grid: { color: gridColor } },
            y: { ticks: { color: textColor, font: { size: 10 } }, grid: { color: gridColor } }
        },
        interaction: { intersect: false, mode: 'index' },
        elements: { point: { radius: 3, hoverRadius: 5 }, line: { tension: 0.35, borderWidth: 2 } }
    };

    // Destroy existing charts
    Object.values(chartInstances).forEach(c => c.destroy());
    chartInstances = {};

    // Conversion Rate
    chartInstances.conv = new Chart(document.getElementById('chartConversion'), {
        type: 'line',
        data: {
            labels: days,
            datasets: [
                { label: 'Control', data: [4.8,4.9,5.0,5.1,4.95,5.0,5.05,5.0,5.02], borderColor: controlColor, backgroundColor: controlColor + '20' },
                { label: 'Treatment', data: [4.9,5.1,5.2,5.3,5.35,5.4,5.45,5.48,5.51], borderColor: treatColor, backgroundColor: treatColor + '20' }
            ]
        },
        options: { ...defaultOpts, scales: { ...defaultOpts.scales, y: { ...defaultOpts.scales.y, ticks: { ...defaultOpts.scales.y.ticks, callback: v => v + '%' } } } }
    });

    // Revenue
    chartInstances.rev = new Chart(document.getElementById('chartRevenue'), {
        type: 'line',
        data: {
            labels: days,
            datasets: [
                { label: 'Control', data: [41.2,42.0,42.5,42.1,42.8,42.3,42.0,42.4,42.3], borderColor: controlColor, backgroundColor: controlColor + '20' },
                { label: 'Treatment', data: [41.5,42.2,42.8,42.6,43.0,42.7,42.5,42.9,42.8], borderColor: treatColor, backgroundColor: treatColor + '20' }
            ]
        },
        options: { ...defaultOpts, scales: { ...defaultOpts.scales, y: { ...defaultOpts.scales.y, ticks: { ...defaultOpts.scales.y.ticks, callback: v => '₹' + v } } } }
    });

    // Applications per user
    chartInstances.apps = new Chart(document.getElementById('chartApplications'), {
        type: 'bar',
        data: {
            labels: days,
            datasets: [
                { label: 'Control', data: [1.75,1.80,1.82,1.78,1.84,1.80,1.83,1.81,1.82], backgroundColor: controlColor + '80', borderRadius: 4, barPercentage: 0.4 },
                { label: 'Treatment', data: [1.82,1.90,1.95,1.98,2.00,2.01,2.02,2.03,2.04], backgroundColor: treatColor + '80', borderRadius: 4, barPercentage: 0.4 }
            ]
        },
        options: defaultOpts
    });

    // Error rate
    chartInstances.err = new Chart(document.getElementById('chartErrors'), {
        type: 'line',
        data: {
            labels: days,
            datasets: [
                { label: 'Control', data: [0.30,0.32,0.29,0.31,0.33,0.30,0.31,0.32,0.31], borderColor: controlColor, backgroundColor: controlColor + '20' },
                { label: 'Treatment', data: [0.31,0.30,0.28,0.30,0.29,0.28,0.30,0.29,0.29], borderColor: treatColor, backgroundColor: treatColor + '20' }
            ]
        },
        options: { ...defaultOpts, scales: { ...defaultOpts.scales, y: { ...defaultOpts.scales.y, ticks: { ...defaultOpts.scales.y.ticks, callback: v => v + '%' } } } }
    });
}
