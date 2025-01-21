#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <tuple>
#include <wx/wx.h>
#include <wx/grid.h>

// Função f(x) e sua derivada
float f(float a3, float a2, float d) {
    return a3 * std::pow(d, 3) - 9.0f * a2 * d + 3.0f;
}

float df(float a3, float a2, float d) {
    return 3.0f * a3 * d * d - 9.0f * a2;
}

// Método de Newton com saída detalhada
std::tuple<float, std::vector<std::tuple<int, float, float, float, float>>> newton(
    float a3, float a2, float x0, float eps, int maxIter) {
    std::vector<std::tuple<int, float, float, float, float>> data;

    for (int i = 0; i < maxIter; ++i) {
        float x1 = x0 - f(a3, a2, x0) / df(a3, a2, x0);
        float fk = std::fabs(f(a3, a2, x1));
        float dx = std::fabs(x1 - x0);
        
        data.emplace_back(i, x1, fk, dx, x0);

        if (fk <= eps || dx <= eps) {
            return {x1, data};
        }

        x0 = x1;
    }

    return {x0, data};
}

std::tuple<float, std::vector<std::tuple<int, float, float, float, float>>> newtonFl(
    float a3, float a2, float x0, float eps, float lambda, int maxIter) {
    std::vector<std::tuple<int, float, float, float, float>> data;

    if(std::fabs(f(a3, a2, x0)) <= eps) {
        return {x0, data};
    }

    float xw = x0;
    for (int i = 0; i < maxIter; ++i) {
        if (std::fabs(df(a3, a2, x0)) >= lambda) {
            xw = x0;
        }

        float fl = std::fabs(df(a3, a2, x0)) > lambda ? df(a3, a2, x0) : df(a3, a2, xw);
        float x1 = x0 - f(a3, a2, x0) / fl;

        float fk = std::fabs(f(a3, a2, x1));

        float dx = std::fabs(x1 - x0);
        
        data.emplace_back(i, x1, fk, dx, x0);

        if (fk <= eps || dx <= eps) {
            return {x1, data};
        }

        x0 = x1;
    }
    return {x0, data};
}

float df2(float a3, float a2, float d, float h = 1e-5f) {
    return (f(a3, a2, d + h) - f(a3, a2, d - h)) / (2 * h);
}

std::tuple<float, std::vector<std::tuple<int, float, float, float, float>>> newtonNumerico(
    float a3, float a2, float x0, float eps, int maxIter) {

    std::vector<std::tuple<int, float, float, float, float>> data;

    if (std::fabs(f(a3, a2, x0)) <= eps) {
        return {x0, data};
    }

    for (int i = 0; i < maxIter; ++i) {
        float x1 = x0 - f(a3, a2, x0) / df2(a3, a2, x0);

        float fk = std::fabs(f(a3, a2, x1));

        float dx = std::fabs(x1 - x0);

        data.emplace_back(i, x1, fk, dx, x0);

        if (fk <= eps || dx <= eps) {
            return {x1, data};
        }

        x0 = x1;
    }

    return {x0, data};
}

// Classe principal da interface wxWidgets
class MyFrame : public wxFrame {
public:
    MyFrame(const wxString& title);

private:
    wxPanel* m_panel;
    wxStaticText* m_title1;
    wxStaticText* m_title2;
    wxStaticText* m_title3;
    wxGrid* m_grid;
    wxGrid* m_grid2;
    wxGrid* m_grid3;
    wxBoxSizer* m_sizer;
    void PopulateGrid(wxGrid* m_grid, const std::vector<std::tuple<int, float, float, float, float>>& data);
};

// Aplicação wxWidgets
class MyApp : public wxApp {
public:
    virtual bool OnInit();
};

wxIMPLEMENT_APP(MyApp);

// Implementação do frame principal
MyFrame::MyFrame(const wxString& title)
    : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(600, 300)) {
    m_panel = new wxPanel(this);

    m_title1 = new wxStaticText(m_panel, wxID_ANY, "Newton Original", wxDefaultPosition, wxDefaultSize, wxALIGN_CENTER);
    m_grid = new wxGrid(m_panel, wxID_ANY);

    m_title2 = new wxStaticText(m_panel, wxID_ANY, "Newton FL", wxDefaultPosition, wxDefaultSize, wxALIGN_CENTER);
    m_grid2 = new wxGrid(m_panel, wxID_ANY);

    m_title3 = new wxStaticText(m_panel, wxID_ANY, "Newton Numerico", wxDefaultPosition, wxDefaultSize, wxALIGN_CENTER);
    m_grid3 = new wxGrid(m_panel, wxID_ANY);

    m_sizer = new wxBoxSizer(wxVERTICAL);

    // Configuração inicial da grid
    m_grid->CreateGrid(0, 5);
    m_grid->SetColLabelValue(0, "Iteracao");
    m_grid->SetColLabelValue(1, "x_(k+1)");
    m_grid->SetColLabelValue(2, "x_k");
    m_grid->SetColLabelValue(3, "f(x_(k+1))");
    m_grid->SetColLabelValue(4, "|x_(k+1) - x_k|");
    m_grid->EnableEditing(false);
    m_grid->SetRowLabelSize(0);
    m_grid->SetColLabelSize(30);


    m_grid->SetColSize(0, 150);
    m_grid->SetColSize(1, 150);
    m_grid->SetColSize(2, 150);
    m_grid->SetColSize(3, 150);
    m_grid->SetColSize(4, 150);

    // Configuração inicial da grid
    m_grid2->CreateGrid(0, 5);
    m_grid2->SetColLabelValue(0, "Iteracao");
    m_grid2->SetColLabelValue(1, "x_(k+1)");
    m_grid2->SetColLabelValue(2, "x_k");
    m_grid2->SetColLabelValue(3, "f(x_(k+1))");
    m_grid2->SetColLabelValue(4, "|x_(k+1) - x_k|");
    m_grid2->EnableEditing(false);
    m_grid2->SetRowLabelSize(0);
    m_grid2->SetColLabelSize(30);


    m_grid2->SetColSize(0, 150);
    m_grid2->SetColSize(1, 150);
    m_grid2->SetColSize(2, 150);
    m_grid2->SetColSize(3, 150);
    m_grid2->SetColSize(4, 150);

    // Configuração inicial da grid
    m_grid3->CreateGrid(0, 5);
    m_grid3->SetColLabelValue(0, "Iteracao");
    m_grid3->SetColLabelValue(1, "x_(k+1)");
    m_grid3->SetColLabelValue(2, "x_k");
    m_grid3->SetColLabelValue(3, "f(x_(k+1))");
    m_grid3->SetColLabelValue(4, "|x_(k+1) - x_k|");
    m_grid3->EnableEditing(false);
    m_grid3->SetRowLabelSize(0);
    m_grid3->SetColLabelSize(30);


    m_grid3->SetColSize(0, 150);
    m_grid3->SetColSize(1, 150);
    m_grid3->SetColSize(2, 150);
    m_grid3->SetColSize(3, 150);
    m_grid3->SetColSize(4, 150);

    // Valores iniciais
    float a3 = 1.0f;
    float a2 = 2.0f;
    float x0 = 1.0f;
    float eps = 1e-6f;
    float lambda = 0.05f;
    int maxIter = 20;

    // Calculando os dados pelo método de Newton
    auto [raiz1, data1] = newton(a3, a2, x0, eps, maxIter);
    auto [raiz2, data2] = newtonFl(a3, a2, x0, eps, lambda, maxIter);
    auto [raiz3, data3] = newtonNumerico(a3, a2, x0, eps, maxIter);

    // Preenchendo a grid com os dados
    PopulateGrid(m_grid, data1);
    PopulateGrid(m_grid2, data2);
    PopulateGrid(m_grid3, data3);

    // Adicionar Tabela 1
    m_sizer->Add(m_title1, 0, wxALIGN_CENTER | wxTOP | wxBOTTOM, 2);
    m_sizer->Add(m_grid, 1, wxEXPAND | wxALL, 2);

    // Adicionar Tabela 2
    m_sizer->Add(m_title2, 0, wxALIGN_CENTER | wxTOP | wxBOTTOM, 2);
    m_sizer->Add(m_grid2, 1, wxEXPAND | wxALL, 2);

    m_sizer->Add(m_title3, 0, wxALIGN_CENTER | wxTOP | wxBOTTOM, 2);
    m_sizer->Add(m_grid3, 1, wxEXPAND | wxALL, 2);

    // Configurar sizer do painel
    m_panel->SetSizer(m_sizer);

    // Ajusta o tamanho da janela ao conteúdo
    m_sizer->Fit(this);
}

// Preencher grid com os dados do método de Newton
void MyFrame::PopulateGrid(wxGrid* m_grid, const std::vector<std::tuple<int, float, float, float, float>>& data) {
    m_grid->AppendRows(data.size());
    int row = 0;
    for (const auto& [k, x1, fk, dx, x0] : data) {
        m_grid->SetCellValue(row, 0, std::to_string(k));
        m_grid->SetCellValue(row, 1, std::to_string(x1));
        m_grid->SetCellValue(row, 2, std::to_string(x0));
        m_grid->SetCellValue(row, 3, std::to_string(fk));
        m_grid->SetCellValue(row, 4, std::to_string(dx));
        ++row;
    }
}

// Inicialização da aplicação wxWidgets
bool MyApp::OnInit() {
    MyFrame* frame = new MyFrame("Newton");
    frame->Show(true);
    return true;
}
