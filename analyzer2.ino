//#define PRINTDATA true
// преобразование Фурье на C (для Raspberry Pi Pico W)
#define FHT_N 1024
#define FFT_SIZE FHT_N

static float sinF[] = {0.0, -1.0, -0.707107, -0.382683, -0.195090, -0.098017, -0.049068, -0.024541, -0.012272, -0.006136};

void FFT(int* AVal, int* FTvl) {
    int i, j, m, Mmax, Istp, count = 0;
    float Tmpr, Tmpi, Tmvl[FFT_SIZE * 2];
    float Wpr, Wr, Wi;

    for (i = 0; i < FFT_SIZE * 2; i += 2) {
        Tmvl[i] = 0;
        Tmvl[i + 1] = AVal[i / 2];
    }

    i = j = 1;
    while (i < FFT_SIZE * 2) {
        if (j > i) {
            Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
            Tmpr = Tmvl[i + 1]; Tmvl[i + 1] = Tmvl[j + 1]; Tmvl[j + 1] = Tmpr;
        }
        i += 2;
        m = FFT_SIZE;
        while ((m >= 2) && (j > m)) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    Mmax = 2;
    while (FFT_SIZE * 2 > Mmax) {
        Wpr = sinF[count + 1] * sinF[count + 1] * 2;
        Istp = Mmax * 2;
        Wr = 1;
        Wi = 0;
        m = 1;

        while (m < Mmax) {
            i = m;
            m += 2;
            Tmpr = Wr; Tmpi = Wi;
            Wr += -Tmpr * Wpr - Tmpi * sinF[count];
            Wi += Tmpr * sinF[count] - Tmpi * Wpr;

            while (i < FFT_SIZE * 2) {
                j = i + Mmax;
                Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j - 1];
                Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j - 1];

                Tmvl[j] = Tmvl[i] - Tmpr;
                Tmvl[j - 1] = Tmvl[i - 1] - Tmpi;
                Tmvl[i] += Tmpr;
                Tmvl[i - 1] += Tmpi;
                i += Istp;
            }
        }
        count++;
        Mmax = Istp;
    }
    
    for (i = 0; i < FFT_SIZE; i++) {
        j = i * 2;
        FTvl[i] = (int)(Tmvl[j] * Tmvl[j] + Tmvl[j + 1] * Tmvl[j + 1]) >> 20;
    }
}

#include <Adafruit_NeoPixel.h>
#define LED_PIN 11

#define SMOOTH_UP 1.0f
#define SMOOTH_DOWN 0.35f
#define BLOCKS 5
#define LEVELS 8
#define COLUMNS (BLOCKS*8)
#define NUM_LEDS (COLUMNS*LEVELS)
#define GAIN 130.f
#define LOW_PASS -100.f
#define Fq 48000.f
#define TS (1.f/Fq*1000000.)

Adafruit_NeoPixel matrix = Adafruit_NeoPixel(NUM_LEDS, LED_PIN, NEO_GRB + NEO_KHZ800);

// Двойная буферизация для параллельной обработки
int raw_buffer1[FFT_SIZE];
int raw_buffer2[FFT_SIZE];
int spectr_buffer1[FFT_SIZE];
int spectr_buffer2[FFT_SIZE];

int* current_raw = raw_buffer1;
int* current_spectr = spectr_buffer1;
int* processing_raw = raw_buffer2;
int* processing_spectr = spectr_buffer2;

uint8_t v_array_max[COLUMNS];
uint8_t v_array[COLUMNS];
float v_array_old[COLUMNS];

#define FALL_DELAY 30
#define FALL_PAUSE 400
unsigned long timeMaxLevel[COLUMNS];
bool fallFlag;
unsigned long fallTimer;

int band_start[COLUMNS];
int band_end[COLUMNS];
float band_compensation[COLUMNS];

// Флаги для синхронизации между ядрами
volatile bool data_ready = false;
volatile bool processing_done = false;
volatile bool new_data_available = false;

void setupFrequencyBands() {
    float min_freq = 90.0f;
    float max_freq = 24000.0f;
    
    float min_log = log10f(min_freq);
    float max_log = log10f(max_freq);
    float log_range = max_log - min_log;
    
    for (int band = 0; band < COLUMNS; band++) {
        float band_pos = (float)band / COLUMNS;
        float next_band_pos = (float)(band + 1) / COLUMNS;
        
        float low_freq = powf(10.0f, min_log + band_pos * log_range);
        float high_freq = powf(10.0f, min_log + next_band_pos * log_range);
        
        band_start[band] = (int)(low_freq * FFT_SIZE / Fq);
        band_end[band] = (int)(high_freq * FFT_SIZE / Fq);
        
        band_start[band] = constrain(band_start[band], 1, FFT_SIZE/2 - 1);
        band_end[band] = constrain(band_end[band], band_start[band] + 1, FFT_SIZE/2);
        
        float center_freq = sqrtf(low_freq * high_freq);
        
        //if (center_freq < 900.0f) {
        //    band_compensation[band] = 1.0f;
        //} else {
            band_compensation[band] = 1.0f + (center_freq+900) / 6000.0f;
        //}
        band_compensation[band] = constrain(band_compensation[band], 1.0f, 1.8f);
    }
}

float mapf(float x, float in_min, float in_max, float out_min, float out_max) {
    return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

void runDisplay() {
    matrix.clear();
    
    for (uint16_t c = 0; c < COLUMNS; ++c) {
        for (uint16_t l = 0; l < LEVELS; ++l) {
            uint16_t idx = NUM_LEDS - (c * LEVELS + LEVELS - l);
            
            if (l <= v_array[c]) {
                matrix.setPixelColor(idx, matrix.Color(0, 1, 1));
            }
            if (l == v_array_max[c]) {
                matrix.setPixelColor(idx, matrix.Color(1, 0, 0));
            }
        }
    }
    matrix.show();
}

void runPrepareArray() {
#ifdef PRINTDATA
    Serial.print("db data :  ");
#endif        

    for (int band = 0; band < COLUMNS; band++) {
        if (band_start[band] >= band_end[band]) {
            v_array[band] = 0;
            continue;
        }
        
        float sum = 0;
        int count = 0;
        
        for (int bin = band_start[band]; bin < band_end[band]; bin++) {
            sum += (float)current_spectr[bin];
            count++;
        }
        
        if (count == 0) {
            v_array[band] = 0;
            continue;
        }
        
        float avg = (sum / count);///4;
        float dB = 20.0f * log10f(avg + 1.0f);
        
#ifdef PRINTDATA
        Serial.print(dB);
        Serial.print(" ");
#endif        
        
        dB *= band_compensation[band];
        
        //float scaled_dB = mapf(dB, 10.0f, 80.0f, LOW_PASS, GAIN);
        //scaled_dB = constrain(scaled_dB, LOW_PASS, GAIN);

        float scaled_dB = dB;//mapf(dB, 10.0f, 80.0f, LOW_PASS, GAIN);
        //scaled_dB = constrain(scaled_dB, LOW_PASS, GAIN);
        
        float smooth_factor = (scaled_dB > v_array_old[band]) ? SMOOTH_UP : SMOOTH_DOWN;
        float smoothed_dB = scaled_dB * smooth_factor + v_array_old[band] * (1.0f - smooth_factor);
        v_array_old[band] = smoothed_dB;
        
        //v_array[band] = (uint8_t)constrain(mapf(smoothed_dB, LOW_PASS, GAIN, 0.0f, (float)LEVELS), 0, LEVELS - 1);
        
        v_array[band] = (uint8_t)/*constrain(*/ mapf(smoothed_dB, 20.f, 90.f, 0.0f, LEVELS-1.f);//, 0, LEVELS-1);
        
        if (v_array[band] >= v_array_max[band]) {
            v_array_max[band] = v_array[band];
            timeMaxLevel[band] = millis();
        }
        
        if (fallFlag && (millis() - timeMaxLevel[band] > FALL_PAUSE) && v_array_max[band] > 0) {
            v_array_max[band]--;
        }
    }

#ifdef PRINTDATA
    Serial.println();
#endif        

    fallFlag = (millis() - fallTimer > FALL_DELAY);
    if (fallFlag) fallTimer = millis();
}

void getData() {
    for (int i = 0; i < FFT_SIZE; i++) {
        current_raw[i] = analogRead(A0) - 2048;
        delayMicroseconds(TS);
    }
}

void setup() {
#ifdef PRINTDATA
    Serial.begin(115200);
    while (!Serial) delay(10);
#endif

    analogReadResolution(12);
    
    matrix.begin();
    matrix.setBrightness(255);
    //matrix.fill(matrix.Color(1, 0, 1));
    //matrix.show();
    //delay(25);
    matrix.clear();
    matrix.show();
    
    setupFrequencyBands();
    
    // Запускаем обработку на втором ядре
    data_ready = true;
}

void loop() {
    // Ядро 0: сбор данных и отображение
    static unsigned long time, tData, tPrepare, tDisplay;
    
    time = micros();
    getData();
    tData = micros() - time;
    
    // Сигнализируем второму ядру, что данные готовы
    data_ready = true;
    
    // Ждем завершения обработки на втором ядре
    while (!processing_done) {
        delayMicroseconds(2);
    }
    
    time = micros();
    runPrepareArray();
    tPrepare = micros() - time;
    
    time = micros();
    runDisplay();
    tDisplay = micros() - time;

#ifdef PRINTDATA
        // Отладочный вывод
        //Serial.print("Times: data: ");
        //Serial.print(tData);
        //Serial.print("us, prepare: ");
        //Serial.print(tPrepare);
        //Serial.print("us, display: ");
        //Serial.print(tDisplay);
        //Serial.print("us, ");
        Serial.print("vs_array: ");
        for (int i = 0; i < COLUMNS; i++) {
            Serial.print(v_array[i]);
            Serial.print(" ");
        }
        Serial.println();
#endif        
    // Сигнализируем, что можно использовать новые данные
    new_data_available = true;
}

void loop1() {
    // Ядро 1: обработка FFT
    static unsigned long tFFT;
    
    while (true) {
        // Ждем новые данные от первого ядра
        while (!data_ready) {
            delayMicroseconds(2);
        }
        
        data_ready = false;
        
        // Выполняем FFT
        unsigned long time = micros();
        FFT(current_raw, current_spectr);
        tFFT = micros() - time;

#ifdef PRINTDATA
        // Отладочный вывод
        Serial.print("raw data: ");
        for (int i = 0; i < COLUMNS; i++) {
            Serial.print(current_raw[i]);
            Serial.print(" ");
        }
        Serial.println();
        Serial.print("spc data: ");
        for (int i = 0; i < COLUMNS; i++) {
            Serial.print(current_spectr[i]);
            Serial.print(" ");
        }
        Serial.println();
#endif        

        // Сигнализируем о завершении обработки
        processing_done = true;
        
        // Ждем, пока первое ядро не обработает результаты
        while (!new_data_available) {
            delayMicroseconds(2);
        }
        
        new_data_available = false;
        processing_done = false;
        
        // Меняем буферы местами для следующей итерации
        int* temp_raw = current_raw;
        int* temp_spectr = current_spectr;
        
        current_raw = processing_raw;
        current_spectr = processing_spectr;
        
        processing_raw = temp_raw;
        processing_spectr = temp_spectr;
#ifdef PRINTDATA
        // Вывод времени FFT для отладки
        //Serial.print("FFT time: ");
        //Serial.print(tFFT);
        //Serial.println("us");
#endif        
    }
}
