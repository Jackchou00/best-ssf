// 全局变量
let LMS = [];
let wavelengths = [];
let max_response = [];
let bestCCM = null;
let bestKMax = 0;
let ssfChart = null;

// 初始化
document.addEventListener('DOMContentLoaded', function() {
    // 加载CSV数据
    loadCSVData();
    
    // 设置按钮事件监听器
    document.getElementById('runButton').addEventListener('click', runIterations);
    document.getElementById('saveButton').addEventListener('click', saveResults);
});

// 加载CSV数据
function loadCSVData() {
    Papa.parse('static/CIE_xyz_2015_10deg.csv', {
        download: true,
        delimiter: ",",
        complete: function(results) {
            // 处理CSV数据
            const data = results.data;
            // 过滤掉空行和转换为数值
            const numericData = data
                .filter(row => row.length > 1)
                .map(row => row.map(val => parseFloat(val)))
                .filter(row => !isNaN(row[0]));
            
            // 只取前391行的数据，并且只取第2列到第4列
            LMS = numericData.slice(0, 391).map(row => row.slice(1, 4));
            
            // 生成波长数据
            wavelengths = Array.from({length: 391}, (_, i) => 390 + i);
            
            // 计算max_response
            max_response = wavelengths.map(w => w / 780);
            
            // 更新状态
            document.getElementById('status').textContent = 'Data loaded. Ready to run iterations.';
            
            // 初始化图表
            initChart();
        },
        error: function(error) {
            console.error('Error loading CSV:', error);
            document.getElementById('status').textContent = 'Error loading data. Check console for details.';
        }
    });
}

// 初始化图表
function initChart() {
    const ctx = document.getElementById('ssfChart').getContext('2d');
    ssfChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: wavelengths,
            datasets: [
                {
                    label: 'Max Response',
                    data: max_response,
                    borderColor: 'rgba(0, 0, 0, 0.5)',
                    borderDash: [5, 5],
                    fill: false
                }
            ]
        },
        options: {
            responsive: true,
            title: {
                display: true,
                text: 'Spectral Sensitivity Function'
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Wavelength (nm)'
                    }
                },
                y: {
                    title: {
                        display: true,
                        text: 'Relative Response'
                    }
                }
            }
        }
    });
}

// 运行迭代
function runIterations() {
    const button = document.getElementById('runButton');
    button.disabled = true;
    document.getElementById('status').textContent = 'Running iterations...';
    
    // 使用setTimeout来避免UI阻塞
    setTimeout(() => {
        generateCCM(10000);
        button.disabled = false;
        document.getElementById('saveButton').disabled = false;
        document.getElementById('status').textContent = 'Iterations complete!';
    }, 100);
}

// 生成CCM矩阵
function generateCCM(maxIter) {
    // 计算V (LMS的列和)
    const V = [0, 0, 0];
    for (let i = 0; i < LMS.length; i++) {
        for (let j = 0; j < 3; j++) {
            V[j] += LMS[i][j];
        }
    }
    
    // 迭代寻找最佳CCM
    for (let iter = 0; iter < maxIter; iter++) {
        // 生成随机CCM矩阵
        const C = [
            [Math.random() * 2 - 1, Math.random() * 2 - 1, Math.random() * 2 - 1],
            [Math.random() * 2 - 1, Math.random() * 2 - 1, Math.random() * 2 - 1],
            [Math.random() * 2 - 1, Math.random() * 2 - 1, Math.random() * 2 - 1]
        ];
        
        // 检查线性无关性 (简化版)
        const det = determinant3x3(C);
        if (Math.abs(det) < 1e-10) continue;
        
        const cList = [];
        let validVectors = true;
        
        for (let j = 0; j < 3; j++) {
            const c = C[j];
            // 计算点积
            let currentDot = 0;
            for (let k = 0; k < 3; k++) {
                currentDot += V[k] * c[k];
            }
            
            // 如果点积接近0，这个向量不可用
            if (Math.abs(currentDot) < 1e-10) {
                validVectors = false;
                break;
            }
            
            // 归一化
            const normalizedC = c.map(val => val / currentDot);
            cList.push(normalizedC);
        }
        
        if (!validVectors) continue;
        
        // 计算k_max
        let kMax = Infinity;
        for (let j = 0; j < 3; j++) {
            const cj = cList[j];
            
            // 计算所有LMS与当前c的点积
            const responses = [];
            for (let i = 0; i < LMS.length; i++) {
                let dotProduct = 0;
                for (let k = 0; k < 3; k++) {
                    dotProduct += LMS[i][k] * cj[k];
                }
                responses.push(dotProduct);
            }
            
            // 只考虑正响应
            const positiveIndices = [];
            for (let i = 0; i < responses.length; i++) {
                if (responses[i] > 0) {
                    positiveIndices.push(i);
                }
            }
            
            if (positiveIndices.length === 0) {
                validVectors = false;
                break;
            }
            
            // 计算比率
            const ratios = positiveIndices.map(i => max_response[i] / responses[i]);
            const minRatio = Math.min(...ratios);
            
            if (minRatio < kMax) {
                kMax = minRatio;
            }
        }
        
        // 跳过无效的情况
        if (!validVectors || kMax === Infinity || kMax <= 0) continue;
        
        // 更新最佳解
        if (kMax > bestKMax) {
            bestKMax = kMax;
            
            // 构建CCM矩阵
            bestCCM = [[], [], []];
            for (let i = 0; i < 3; i++) {
                for (let j = 0; j < 3; j++) {
                    bestCCM[i][j] = cList[j][i] * kMax;
                }
            }
        }
    }
    
    // 验证结果并更新UI
    validateAndUpdateUI();
}

// 验证结果并更新UI
function validateAndUpdateUI() {
    if (!bestCCM) {
        document.getElementById('ccmMatrix').textContent = 'Failed to generate valid CCM matrix';
        return;
    }
    
    // 显示CCM矩阵
    document.getElementById('ccmMatrix').textContent = 
        bestCCM.map(row => row.map(val => val.toFixed(6)).join('\t')).join('\n');
    
    // 计算行列式
    const det = determinant3x3(bestCCM);
    document.getElementById('determinant').textContent = det.toFixed(6);
    
    // 计算SSF
    const SSF = multiplyMatrices(LMS, bestCCM);
    
    // 检查约束条件
    let constraintMet = true;
    for (let i = 0; i < max_response.length; i++) {
        if (SSF[i][0] > max_response[i] || 
            SSF[i][1] > max_response[i] || 
            SSF[i][2] > max_response[i]) {
            constraintMet = false;
            break;
        }
    }
    document.getElementById('constraintsMet').textContent = constraintMet ? 'Yes' : 'No';
    
    // 计算V (LMS的列和)
    const V = [0, 0, 0];
    for (let i = 0; i < LMS.length; i++) {
        for (let j = 0; j < 3; j++) {
            V[j] += LMS[i][j];
        }
    }
    
    // 计算RGB
    const RGB = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            RGB[i] += V[j] * bestCCM[j][i];
        }
    }
    
    document.getElementById('rgbValues').textContent = 
        RGB.map(val => val.toFixed(6)).join(', ');
    
    // 检查RGB是否相等
    const rgbEqual = Math.abs(RGB[0] - RGB[1]) < 1e-6 && Math.abs(RGB[1] - RGB[2]) < 1e-6;
    document.getElementById('rgbEqual').textContent = rgbEqual ? 'Yes' : 'No';
    
    // 检查CCM是否可逆
    document.getElementById('ccmInvertible').textContent = Math.abs(det) > 1e-10 ? 'Yes' : 'No';
    
    // 更新图表
    updateChart(SSF);
}

// 更新图表
function updateChart(SSF) {
    // 提取SSF的三个通道
    const channel1 = SSF.map(row => row[0]);
    const channel2 = SSF.map(row => row[1]);
    const channel3 = SSF.map(row => row[2]);
    
    // 更新图表数据
    if (ssfChart.data.datasets.length === 1) {
        // 添加新的数据集
        ssfChart.data.datasets.push(
            {
                label: 'Channel 1',
                data: channel1,
                borderColor: 'rgba(255, 0, 0, 0.8)',
                fill: false
            },
            {
                label: 'Channel 2',
                data: channel2,
                borderColor: 'rgba(0, 255, 0, 0.8)',
                fill: false
            },
            {
                label: 'Channel 3',
                data: channel3,
                borderColor: 'rgba(0, 0, 255, 0.8)',
                fill: false
            }
        );
    } else {
        // 更新现有数据集
        ssfChart.data.datasets[1].data = channel1;
        ssfChart.data.datasets[2].data = channel2;
        ssfChart.data.datasets[3].data = channel3;
    }
    
    ssfChart.update();
}

// 保存结果
function saveResults() {
    if (!bestCCM) {
        alert('No valid CCM matrix to save');
        return;
    }
    
    // 创建要保存的文本内容
    let content = 'CCM Matrix:\n';
    content += bestCCM.map(row => row.map(val => val.toFixed(6)).join(',')).join('\n');
    content += '\n\nValidation Results:\n';
    content += `Determinant: ${document.getElementById('determinant').textContent}\n`;
    content += `Constraints Met: ${document.getElementById('constraintsMet').textContent}\n`;
    content += `RGB Values: ${document.getElementById('rgbValues').textContent}\n`;
    content += `RGB Equal: ${document.getElementById('rgbEqual').textContent}\n`;
    content += `CCM Invertible: ${document.getElementById('ccmInvertible').textContent}\n`;
    
    // 创建Blob对象
    const blob = new Blob([content], { type: 'text/plain' });
    
    // 创建下载链接
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'ccm_results.txt';
    
    // 触发下载
    document.body.appendChild(a);
    a.click();
    
    // 清理
    document.body.removeChild(a);
    URL.revokeObjectURL(a.href);
}

// 辅助函数：3x3矩阵行列式
function determinant3x3(matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

// 辅助函数：矩阵乘法
function multiplyMatrices(A, B) {
    const result = [];
    
    for (let i = 0; i < A.length; i++) {
        result[i] = [];
        for (let j = 0; j < B[0].length; j++) {
            let sum = 0;
            for (let k = 0; k < A[0].length; k++) {
                sum += A[i][k] * B[k][j];
            }
            result[i][j] = sum;
        }
    }
    
    return result;
}
