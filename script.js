/* ============================================
   IIMJobs Experiment Hub — Application Logic v2.1
   ============================================ */

document.addEventListener('DOMContentLoaded', () => {
    initNavigation();
    initThemeToggle();
    initMobileMenu();
    buildTermsGlossary();
});

/* ========================
   NAVIGATION
   ======================== */
function initNavigation() {
    const navItems = document.querySelectorAll('.nav-item');
    navItems.forEach(item => {
        item.addEventListener('click', (e) => {
            e.preventDefault();
            const sectionId = item.getAttribute('data-section');
            switchSection(sectionId);
            document.getElementById('sidebar').classList.remove('open');
        });
    });
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
}

/* ========================
   THEME TOGGLE
   ======================== */
function initThemeToggle() {
    const stored = localStorage.getItem('theme') || 'dark';
    document.documentElement.setAttribute('data-theme', stored);
    document.getElementById('themeToggle').addEventListener('click', toggleTheme);
    const mobileToggle = document.getElementById('themeToggleMobile');
    if (mobileToggle) mobileToggle.addEventListener('click', toggleTheme);
    updateMobileThemeIcon(stored);
}

function toggleTheme() {
    const current = document.documentElement.getAttribute('data-theme');
    const next = current === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    localStorage.setItem('theme', next);
    updateMobileThemeIcon(next);
}

function updateMobileThemeIcon(theme) {
    const btn = document.getElementById('themeToggleMobile');
    if (btn) btn.textContent = theme === 'dark' ? '🌙' : '☀️';
}

/* ========================
   MOBILE MENU
   ======================== */
function initMobileMenu() {
    document.getElementById('hamburger').addEventListener('click', () => {
        document.getElementById('sidebar').classList.toggle('open');
    });
    document.getElementById('mainContent').addEventListener('click', () => {
        document.getElementById('sidebar').classList.remove('open');
    });
}

/* ========================
   STATISTICAL HELPERS
   ======================== */

// Standard normal CDF (Abramowitz & Stegun)
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
    const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
    const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
    const pLow = 0.02425, pHigh = 1 - pLow;
    let q, r;
    if (p < pLow) {
        q = Math.sqrt(-2 * Math.log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    } else if (p <= pHigh) {
        q = p - 0.5; r = q * q;
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    } else {
        q = Math.sqrt(-2 * Math.log(1 - p));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
}

// Log-gamma (Lanczos)
function lnGamma(z) {
    const g = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.001208650973866179, -0.000005395239384953];
    let x = z, y = z, tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += g[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
}

// Regularized lower incomplete gamma (series expansion)
function gammaPSeries(a, x) {
    if (x <= 0) return 0;
    let sum = 1 / a, term = 1 / a;
    for (let n = 1; n < 500; n++) {
        term *= x / (a + n);
        sum += term;
        if (Math.abs(term) < 1e-14) break;
    }
    return Math.min(1, Math.max(0, sum * Math.exp(-x + a * Math.log(x) - lnGamma(a))));
}

// Chi-square CDF
function chiSqCDF(x, df) {
    if (x <= 0) return 0;
    return gammaPSeries(df / 2, x / 2);
}

// Input validation helper
function validateField(id) {
    const el = document.getElementById(id);
    const val = parseFloat(el.value);
    if (isNaN(val) || el.value.trim() === '') {
        el.classList.add('error');
        return null;
    }
    el.classList.remove('error');
    return val;
}

function clearErrors(...ids) {
    ids.forEach(id => document.getElementById(id).classList.remove('error'));
}

function formatNum(n) {
    return n >= 1000000 ? (n / 1000000).toFixed(1) + 'M' : n.toLocaleString('en-US');
}

// Advanced toggle
function toggleAdvanced(btn) {
    const content = btn.nextElementSibling;
    content.classList.toggle('show');
    const icon = btn.querySelector('.adv-icon');
    icon.textContent = content.classList.contains('show') ? '−' : '+';
}

/* ========================
   SAMPLE SIZE CALCULATOR
   ======================== */
let ssGroupCount = 2;

function selectSSGroup(n) {
    ssGroupCount = n;
    document.querySelectorAll('#ssGroupSelect .radio-card').forEach((c, i) => {
        c.classList.toggle('selected', (i === 0 && n === 2) || (i === 1 && n === 3));
    });
    document.getElementById('ss2group').classList.toggle('hidden', n !== 2);
    document.getElementById('ss3group').classList.toggle('hidden', n !== 3);
    document.getElementById('ssResults').classList.add('hidden');
}

function calcSampleSizePair(p1, p2, alpha, power) {
    const zAlpha = normalInv(1 - alpha / 2);
    const zBeta = normalInv(power);
    const pBar = (p1 + p2) / 2;
    const numerator = Math.pow(
        zAlpha * Math.sqrt(2 * pBar * (1 - pBar)) +
        zBeta * Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2)),
        2
    );
    return Math.ceil(numerator / Math.pow(p1 - p2, 2));
}

function calculateSampleSize() {
    const confidence = validateField('ssConfidence');
    const power = validateField('ssPower');
    const traffic = validateField('ssDailyTraffic');
    if (confidence === null || power === null || traffic === null) return;

    const alpha = 1 - confidence / 100;
    const pow = power / 100;
    let perGroup, numGroups;

    if (ssGroupCount === 2) {
        const baseline = validateField('ssBaseline');
        const treatment = validateField('ssTreatment');
        if (baseline === null || treatment === null) return;
        if (baseline === treatment) {
            alert('Baseline and treatment rates must be different.');
            return;
        }
        perGroup = calcSampleSizePair(baseline / 100, treatment / 100, alpha, pow);
        numGroups = 2;
    } else {
        const ctrl = validateField('ssCtrl3');
        const t1 = validateField('ssT1');
        const t2 = validateField('ssT2');
        if (ctrl === null || t1 === null || t2 === null) return;
        // Bonferroni correction: 3 pairwise comparisons
        const adjAlpha = alpha / 3;
        const n1 = calcSampleSizePair(ctrl / 100, t1 / 100, adjAlpha, pow);
        const n2 = calcSampleSizePair(ctrl / 100, t2 / 100, adjAlpha, pow);
        const n3 = calcSampleSizePair(t1 / 100, t2 / 100, adjAlpha, pow);
        perGroup = Math.max(n1, n2, n3);
        numGroups = 3;
    }

    const total = perGroup * numGroups;
    const days = Math.ceil(perGroup / traffic);

    document.getElementById('ssResPerGroup').textContent = formatNum(perGroup);
    document.getElementById('ssResTotal').textContent = formatNum(total);
    document.getElementById('ssResDays').textContent = days + (days === 1 ? ' day' : ' days');
    document.getElementById('ssResults').classList.remove('hidden');
}

/* ========================
   SRM CHECKER
   ======================== */
let srmGroupCount = 2;

function selectSRMGroup(n) {
    srmGroupCount = n;
    document.querySelectorAll('#srmGroupSelect .radio-card').forEach((c, i) => {
        c.classList.toggle('selected', (i === 0 && n === 2) || (i === 1 && n === 3));
    });
    document.getElementById('srm2group').classList.toggle('hidden', n !== 2);
    document.getElementById('srm3group').classList.toggle('hidden', n !== 3);
    document.getElementById('srmResults').classList.add('hidden');
}

function checkSRM() {
    let observed = [], expected = [];

    if (srmGroupCount === 2) {
        const expCtrl = validateField('srmExpCtrl');
        const expTreat = validateField('srmExpTreat');
        const obsCtrl = validateField('srmObsCtrl');
        const obsTreat = validateField('srmObsTreat');
        if ([expCtrl, expTreat, obsCtrl, obsTreat].includes(null)) return;
        const total = obsCtrl + obsTreat;
        const expSum = expCtrl + expTreat;
        expected = [total * (expCtrl / expSum), total * (expTreat / expSum)];
        observed = [obsCtrl, obsTreat];
    } else {
        const expCtrl = validateField('srmExpCtrl3');
        const expT1 = validateField('srmExpT13');
        const expT2 = validateField('srmExpT23');
        const obsCtrl = validateField('srmObsCtrl3');
        const obsT1 = validateField('srmObsT13');
        const obsT2 = validateField('srmObsT23');
        if ([expCtrl, expT1, expT2, obsCtrl, obsT1, obsT2].includes(null)) return;
        const total = obsCtrl + obsT1 + obsT2;
        const expSum = expCtrl + expT1 + expT2;
        expected = [
            total * (expCtrl / expSum),
            total * (expT1 / expSum),
            total * (expT2 / expSum)
        ];
        observed = [obsCtrl, obsT1, obsT2];
    }

    // Chi-square test
    let chiSq = 0;
    for (let i = 0; i < observed.length; i++) {
        chiSq += Math.pow(observed[i] - expected[i], 2) / expected[i];
    }
    const df = observed.length - 1;
    const pValue = 1 - gammaPSeries(df / 2, chiSq / 2);

    // Display badge
    const badge = document.getElementById('srmBadge');
    if (pValue >= 0.05) {
        badge.innerHTML = '<div class="status-badge healthy">🟢 Traffic split looks healthy. No mismatch detected.</div>';
    } else {
        badge.innerHTML = '<div class="status-badge warning">🔴 Traffic split mismatch detected. Experiment may be invalid.</div>';
    }

    // Advanced stats
    document.getElementById('srmAdvanced').innerHTML =
        `<strong>Chi-square statistic:</strong> ${chiSq.toFixed(4)}<br>` +
        `<strong>Degrees of freedom:</strong> ${df}<br>` +
        `<strong>p-value:</strong> ${pValue < 0.0001 ? pValue.toExponential(4) : pValue.toFixed(6)}<br><br>` +
        `If the p-value is below 0.05, there is a statistically significant mismatch between the expected and observed traffic split.`;

    document.getElementById('srmResults').classList.remove('hidden');
}

/* ========================
   LIFT CALCULATOR
   ======================== */
function calculateLift() {
    const ctrlUsers = validateField('liftCtrlUsers');
    const ctrlConv = validateField('liftCtrlConv');
    const treatUsers = validateField('liftTreatUsers');
    const treatConv = validateField('liftTreatConv');
    const confidence = validateField('liftConfidence');
    if ([ctrlUsers, ctrlConv, treatUsers, treatConv, confidence].includes(null)) return;
    if (ctrlUsers <= 0 || treatUsers <= 0) return;

    const pC = ctrlConv / ctrlUsers;
    const pT = treatConv / treatUsers;
    const absLift = pT - pC;
    const relLift = pC > 0 ? absLift / pC : 0;

    // Unpooled SE for CI
    const se = Math.sqrt(pC * (1 - pC) / ctrlUsers + pT * (1 - pT) / treatUsers);
    const z = se > 0 ? absLift / se : 0;
    const pValue = 2 * (1 - normalCDF(Math.abs(z)));

    const zCI = normalInv(1 - (1 - confidence / 100) / 2);
    const ciLow = absLift - zCI * se;
    const ciHigh = absLift + zCI * se;

    // Display result cards
    document.getElementById('liftCtrlRate').textContent = (pC * 100).toFixed(2) + '%';
    document.getElementById('liftTreatRate').textContent = (pT * 100).toFixed(2) + '%';

    const absEl = document.getElementById('liftAbsLift');
    absEl.textContent = (absLift >= 0 ? '+' : '') + (absLift * 100).toFixed(2) + 'pp';
    absEl.style.color = absLift >= 0 ? 'var(--success)' : 'var(--danger)';

    const relEl = document.getElementById('liftRelLift');
    relEl.textContent = (relLift >= 0 ? '+' : '') + (relLift * 100).toFixed(2) + '%';
    relEl.style.color = relLift >= 0 ? 'var(--success)' : 'var(--danger)';

    // CI text
    const confPct = confidence.toFixed(0);
    document.getElementById('liftCIText').textContent =
        `We are ${confPct}% confident that the true lift lies between:`;
    document.getElementById('liftCIRange').textContent =
        `${(ciLow * 100).toFixed(2)}%  and  ${(ciHigh * 100 >= 0 ? '+' : '')}${(ciHigh * 100).toFixed(2)}%`;

    // CI bar visualization
    const maxAbs = Math.max(Math.abs(ciLow * 100), Math.abs(ciHigh * 100), 0.5);
    const scale = 45 / maxAbs;

    const bar = document.getElementById('liftCIBar');
    const point = document.getElementById('liftCIPoint');
    const leftPct = 50 + ciLow * 100 * scale;
    const rightPct = 50 + ciHigh * 100 * scale;
    bar.style.left = leftPct + '%';
    bar.style.width = Math.max(0, rightPct - leftPct) + '%';

    let barColor;
    if (ciLow > 0) barColor = 'var(--success)';
    else if (ciHigh < 0) barColor = 'var(--danger)';
    else barColor = 'var(--accent)';
    bar.style.background = barColor;
    point.style.background = barColor;
    point.style.left = (50 + absLift * 100 * scale) + '%';

    // Significance badge
    const sigEl = document.getElementById('liftSignificance');
    if (ciLow > 0) {
        sigEl.innerHTML = '<div class="status-badge healthy">🟢 Statistically significant positive result.</div>';
    } else if (ciHigh < 0) {
        sigEl.innerHTML = '<div class="status-badge warning">🔴 Statistically significant negative impact.</div>';
    } else {
        sigEl.innerHTML = '<div class="status-badge neutral">⚠️ Result is NOT statistically significant. The confidence interval crosses zero.</div>';
    }

    // Advanced stats
    document.getElementById('liftAdvanced').innerHTML =
        `<strong>Standard Error:</strong> ${se.toFixed(6)}<br>` +
        `<strong>Z-score:</strong> ${z.toFixed(4)}<br>` +
        `<strong>p-value:</strong> ${pValue < 0.0001 ? pValue.toExponential(4) : pValue.toFixed(6)}<br>` +
        `<strong>${confPct}% CI (Absolute):</strong> [${(ciLow * 100).toFixed(4)}%, ${(ciHigh * 100).toFixed(4)}%]`;

    document.getElementById('liftResults').classList.remove('hidden');
}

/* ========================
   GLOSSARY (Understand the Terms)
   ======================== */
const termsData = {
    Framework: [
        { term: 'Primary Metric', desc: 'The direct measurement you expect the experiment to change, such as click-through rate or apply rate.', advanced: 'Also called the "success metric." Statistical tests are applied to this metric to determine if the treatment had a real effect.' },
        { term: 'Outcome Metric', desc: 'The broader business goal you hope the experiment will improve — such as revenue, total applications, or engagement.', advanced: 'Outcome metrics help determine if primary metric improvements translate into real business value. They are monitored but not used to determine statistical significance.' },
        { term: 'Control', desc: 'The group of users who see the existing (unchanged) experience. This is your baseline for comparison.', advanced: 'The control group provides a counterfactual — what would have happened without the change.' },
        { term: 'Treatment', desc: 'The group of users who see the new (changed) experience being tested.', advanced: 'In a multi-variant test, you may have more than one treatment group (e.g., T1, T2), each testing a different change.' },
        { term: 'Randomization', desc: 'Users are randomly assigned to groups so the comparison is fair and unbiased.', advanced: 'Ensures differences between groups are due to the treatment, not pre-existing user differences. Common methods include user-ID-based hashing.' },
        { term: 'ITT (Intent to Treat)', desc: 'Analyzing all users based on their assigned group, even if they didn\'t fully interact with the experiment.', advanced: 'ITT analysis avoids selection bias. If a user was assigned to treatment but never saw the change, they\'re still counted in treatment.' },
    ],
    Sample: [
        { term: 'Baseline', desc: 'The current conversion rate before you start the experiment.', advanced: 'Typically measured from historical data over a recent, representative time period (e.g., past 2–4 weeks).' },
        { term: 'MDE (Minimum Detectable Effect)', desc: 'The smallest improvement you want to be able to detect.', advanced: 'Smaller MDE requires a larger sample size. Choose based on what\'s practically meaningful to the business.' },
        { term: 'Confidence Level', desc: 'How sure you want to be that the result is real and not random noise. Usually 95%.', advanced: '95% confidence means a 5% chance (α = 0.05) of concluding there\'s an effect when there isn\'t (false positive / Type I error).' },
        { term: 'Power', desc: 'The probability of detecting a real improvement if one truly exists. Usually 80%.', advanced: '80% power means a 20% chance of missing a real effect (β = 0.20, a false negative / Type II error). Higher power requires more users.' },
        { term: 'Sample Size', desc: 'The total number of users needed in each group to produce reliable results.', advanced: 'Depends on baseline rate, MDE, confidence level, and power. Undersized experiments risk inconclusive results.' },
    ],
    SRM: [
        { term: 'Expected Split', desc: 'The planned traffic allocation between groups (e.g., 50/50).', advanced: 'Configured in your experimentation platform. Determines what proportion of eligible users should land in each group.' },
        { term: 'Observed Split', desc: 'The actual number of users that ended up in each group.', advanced: 'Small deviations are normal. The SRM test determines if the deviation is large enough to be concerning.' },
        { term: 'Chi-square', desc: 'A statistical test comparing what we expected to see versus what actually happened.', advanced: 'Calculated as Σ((Observed − Expected)² / Expected). Larger values indicate greater discrepancy.' },
        { term: 'p-value', desc: 'A number showing how likely the difference is due to random chance. Below 0.05 means something is probably wrong.', advanced: 'In SRM checks, p < 0.05 means the split is significantly different from planned, suggesting a bug.' },
        { term: 'SRM (Sample Ratio Mismatch)', desc: 'When actual traffic split doesn\'t match the plan. If this happens, results may be invalid.', advanced: 'Common causes: redirect issues, bot filters affecting groups differently, platform bugs in randomization.' },
    ],
    Lift: [
        { term: 'Absolute Lift', desc: 'The simple difference in rates between treatment and control (in percentage points).', advanced: 'Example: Control = 5%, Treatment = 6% → Absolute lift = 1pp. Shows raw magnitude of change.' },
        { term: 'Relative Lift', desc: 'How much better (or worse) the treatment performed compared to control, as a percentage.', advanced: 'Example: Control = 5%, Treatment = 6% → Relative lift = 20%. Useful for comparing across different baselines.' },
        { term: 'Confidence Interval', desc: 'A range of values that likely contains the true effect of the experiment.', advanced: 'A 95% CI means if we ran this experiment many times, 95% of calculated intervals would contain the true effect.' },
        { term: 'Statistical Significance', desc: 'Whether the result is likely real, or could be due to random chance.', advanced: 'Significant when the p-value is below the threshold (typically 0.05) or when the confidence interval doesn\'t contain zero.' },
        { term: 'Standard Error', desc: 'A measure of how much the estimated lift might vary due to random sampling.', advanced: 'SE decreases as sample size increases. Used to calculate CI: estimate ± Z × SE.' },
    ]
};

const termContainerMap = {
    Framework: 'termsFramework',
    Sample: 'termsSample',
    SRM: 'termsSRM',
    Lift: 'termsLift'
};

function buildTermsGlossary() {
    for (const [group, items] of Object.entries(termsData)) {
        const container = document.getElementById(termContainerMap[group]);
        if (!container) continue;
        items.forEach(item => {
            const div = document.createElement('div');
            div.className = 'accordion';
            div.innerHTML = `
                <button class="accordion-header" onclick="this.parentElement.classList.toggle('open')">
                    ${item.term}
                    <span class="accordion-chevron">▼</span>
                </button>
                <div class="accordion-body">
                    <div class="accordion-content">
                        ${item.desc}
                        ${item.advanced ? `
                            <button class="advanced-toggle" onclick="toggleAdvanced(this)">
                                <span class="adv-icon">+</span> Advanced Explanation
                            </button>
                            <div class="advanced-content">${item.advanced}</div>
                        ` : ''}
                    </div>
                </div>
            `;
            container.appendChild(div);
        });
    }
}
