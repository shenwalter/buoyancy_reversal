// Created by Petr Karnakov on 13.06.2022
// Copyright 2022 ETH Zurich

#include <iostream>
#include <limits>

#include <dump/vtk.h>
#include <func/init_u.h>
#include <func/init_vel.h>
#include <kernel/hydro.h>
#include <parse/config.h>
#include <util/format.h>
#include <util/posthook.h>

template <class M>
struct HookContext {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  FieldCell<Scal> fc_chi;
  FieldCell<Scal> fc_buo;
};

template <class M>
void InitTracerFields(
    const std::vector<FieldCell<typename M::Scal>*>& vfcu,
    const std::vector<std::string>& names, const Vars& var, M& m) {
  auto sem = m.GetSem(__func__);
  struct {
    Vars vart;
  } * ctx(sem);
  auto& t = *ctx;
  for (size_t l = 0; l < vfcu.size(); ++l) {
    auto& fcu = *vfcu[l];
    const std::string prefix = "init_" + names[l];
    if (sem("var" + names[l])) {
      fcu.Reinit(m);
      t.vart.String.Set("init_vf", "list");
      t.vart.String.Set("list_path", var.String[prefix + "_list_path"]);
      t.vart.Int.Set("dim", var.Int["dim"]);
      t.vart.Int.Set("list_ls", 3);
    }
    if (sem.Nested("field" + names[l])) {
      InitVf(fcu, t.vart, m, true);
    }
    if (sem("factor" + names[l])) {
      auto k = var.Double(prefix + "_factor", 1);
      for (auto c : m.AllCells()) {
        fcu[c] *= k;
      }
    }
  }
}

template <class Scal>
Scal EvalPiecewise(
    Scal x, const std::vector<Scal>& xx, const std::vector<Scal>& values) {
  fassert_equal(values.size(), xx.size());
  if (xx.size() == 0) {
    return GetNan<Scal>();
  }
  size_t i = 0;
  while (i < xx.size() && xx[i] <= x) {
    ++i;
  }
  if (i == 0) { // x < xx[0]
    return values.front();
  }
  if (i < xx.size()) { // xx[i - 1] <= x < xx[i]
    const Scal x0 = xx[i - 1];
    const Scal x1 = xx[i];
    const Scal y0 = values[i - 1];
    const Scal y1 = values[i];
    return x0 < x1 ? y0 + (y1 - y0) * (x - x0) / (x1 - x0) : y0;
  } else {
    return values.back();
  }
}

// Computes stream function for cloud forcing (Brenguier, Grabowski 1993).
// from Grabowski 1993 Cumulus entrainment ...
template <class M>
void CalcPsiCloud(
    FieldCell<typename M::Scal>& fc_psi, typename M::Vect center,
    typename M::Vect unit, int niter, M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const Scal pi = M_PI;
  fc_psi.Reinit(m, 0);
  for (auto c : m.AllCells()) {
    const Vect xc = (m.GetCenter(c) - center) / unit;
    const Scal x = std::abs(xc[0]);
    const Scal z = xc[1];
    auto& psi = fc_psi[c];
    if (z > 0 && z < 1) {
      using std::exp;
      using std::pow;
      using std::sin;
      psi = z * (z - 1);
      for (int i = 0; i < niter; ++i) {
        const Scal q = 2 * i + 1;
        const Scal dpsi =
            pow(2 / pi, 3) * sin(q * pi * z) / pow(q, 3) * exp(-q * pi * x);
        psi += dpsi;
      }
      if (xc[0] < 0) {
        psi = -psi;
      }
    }
  }
}

template <class M>
void InitVelHook(
    FieldCell<typename M::Vect>& fcv, Hydro<M>* hydro, const Vars& var,
    const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const Vect unit(var.Vect["psi_unit"]);
  const Scal unit_u = var.Double["psi_unit_u"];
  FieldCell<Scal> fc_psi;
  CalcPsiCloud(
      fc_psi, Vect(m.GetGlobalLength()[0] * 0.5, 0), unit, var.Int["psi_niter"],
      m);
  const FieldFace<Scal> ff_psi = UEmbed<M>::Interpolate(
      fc_psi, GetBCondZeroGrad<Scal>(hydro->mebc_fluid_), m);
  const FieldCell<Vect> fc_gpsi = UEmbed<M>::Gradient(ff_psi, m);
  for (auto c : m.SuCells()) {
    fcv[c][0] = fc_gpsi[c][1] * unit_u * unit[1];
    fcv[c][1] = -fc_gpsi[c][0] * unit_u * unit[1];
  }
}

template <class M>
void InitHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  auto sem = m.GetSem();
  if (sem()) {
    hydro->hook_context_ =
        std::make_unique<Holder<HookContext<M>>>(new HookContext<M>());
    auto* hctx =
        dynamic_cast<Holder<HookContext<M>>*>(hydro->hook_context_.get())
            ->Get();
    auto& ht = *hctx;
    ht.fc_chi.Reinit(m, 0);
    ht.fc_buo.Reinit(m, 0);
  }
  if (sem.Nested()) {
    auto* hctx =
        dynamic_cast<Holder<HookContext<M>>*>(hydro->hook_context_.get())
            ->Get();
    auto& ht = *hctx;
    InitTracerFields<M>(
        {
            &ht.fc_chi,
        },
        {
            "chi",
        },
        hydro->var, m);
  }
}

template <class M>
void StepHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  auto& m = hydro->m;
  auto sem = m.GetSem();
  auto& var = hydro->var;
  auto* hctx =
      dynamic_cast<Holder<HookContext<M>>*>(hydro->hook_context_.get())->Get();
  auto& ht = *hctx;
  if (sem.Nested() && hydro->dumper_.Try(hydro->st_.t, hydro->st_.dt)) {
    m.Dump(&ht.fc_chi, "chi");
    m.Dump(&ht.fc_buo, "buo");
  }
  if (sem()) {
    const auto& eb = m;
    using UEB = UEmbed<M>;
    const std::array<FieldCell<Scal>*, 1> fields{
        &ht.fc_chi,
    };
    const Scal diffusion = var.Double["buo_diffusion"];
    for (size_t i = 0; i < fields.size(); ++i) {
      auto pfcu = fields[i];
      auto& fcu = *pfcu;
      const Scal dt = hydro->fs_->GetTimeStep();
      const auto mebc = GetBCondZeroGrad<Scal>(hydro->mebc_fluid_);
      const auto ffg = UEB::Gradient(fcu, mebc, eb);
      const auto fcg = UEB::AverageGradient(ffg, eb);
      const auto& fev = hydro->fs_->GetVolumeFlux();
      const auto feu = UEmbed<M>::InterpolateUpwind(
          fcu, mebc, ConvSc::superbee, fcg, fev, eb);

      FieldEmbed<Scal> fe_flux(m, 0);
      eb.LoopFaces([&](auto cf) {
        // Advection.
        fe_flux[cf] = -feu[cf] * fev[cf];
        // Diffusion.
        fe_flux[cf] += diffusion * ffg[cf] * eb.GetArea(cf);
      });
      for (auto c : eb.CellsM()) {
        Scal sum = 0;
        eb.LoopNci(c, [&](auto q) {
          const auto cf = eb.GetFace(c, q);
          sum += fe_flux[cf] * eb.GetOutwardFactor(c, q);
        });
        fcu[c] += dt * sum / eb.GetVolume(c);
      }
      m.Comm(&fcu);
    }
  }
}

template <class Scal, class T>
T Linear(Scal x, Scal xa, T ua, Scal xb, T ub) {
  if (xa == xb) {
    return (ua + ub) * 0.5;
  }
  return (ua * (xb - x) + ub * (x - xa)) / (xb - xa);
}

// From Fig. 1 (curve 4).
template <class Scal>
Scal BuoF(Scal chi, Scal bc, Scal chis, Scal bs) {
  return chi <= 0     ? bc
         : chi >= 1   ? 0
         : chi < chis ? Linear<Scal, Scal>(chi, 0, bc, chis, bs)
                      : Linear<Scal, Scal>(chi, chis, bs, 1, 0);
}

template <class M>
void PreStepHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto& m = hydro->m;
  auto sem = m.GetSem();
  auto* hctx =
      dynamic_cast<Holder<HookContext<M>>*>(hydro->hook_context_.get())->Get();
  auto& var = hydro->var;
  auto& ht = *hctx;
  if (sem()) {
    const Scal factor = var.Double["buo_factor"];
    const auto bc_y = var.Vect["buo_bc_y"];
    const auto bc_values = var.Vect["buo_bc_values"];
    const auto chis_y = var.Vect["buo_chis_y"];
    const auto chis_values = var.Vect["buo_chis_values"];
    const auto bs_y = var.Vect["buo_bs_y"];
    const auto bs_values = var.Vect["buo_bs_values"];
    auto& fc_force = hydro->fc_force_;
    for (auto c : m.Cells()) {
      const Scal y = m.GetCenter(c)[1];
      const Scal chi = Reconst<Scal>::Clip(ht.fc_chi[c], 0, 1);
      const Scal bc = EvalPiecewise(y, bc_y, bc_values);
      const Scal chis = EvalPiecewise(y, chis_y, chis_values);
      const Scal bs = EvalPiecewise(y, bs_y, bs_values);
      const Scal b = BuoF<Scal>(chi, bc, chis, bs);
      ht.fc_buo[c] = b;
      fc_force[c] += Vect::GetUnit(1) * (b * factor);
    }
  }
}

using M = MeshCartesian<double, 2>;

template void InitHook(Hydro<M>*);
template void InitVelHook(
    FieldCell<typename M::Vect>&, Hydro<M>*, const Vars&, const M&);
template void StepHook(Hydro<M>*);
template void PreStepHook(Hydro<M>*);
